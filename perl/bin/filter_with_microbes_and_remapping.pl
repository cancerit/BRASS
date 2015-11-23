#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;

use DateTime;
use warnings FATAL => 'all';
use Bio::DB::Sam;
use File::Temp qw(tempfile);
use IO::Handle;
use List::Util qw(max sum);
use Getopt::Long;
use Const::Fast qw(const);
use Data::Dumper;
use Capture::Tiny qw(capture);
use File::Path qw(make_path remove_tree);

use Sanger::CGP::Brass::Implement;

const my $DNA_MAT => join("\n",
                      join("\t", q{}, qw(A C G T)),
                      join("\t", q{A}, qw(2 -2 -2 -2)),
                      join("\t", q{C}, qw(-2 2 -2 -2)),
                      join("\t", q{G}, qw(-2 -2 2 -2)),
                      join("\t", q{T}, qw(-2 -2 -2 2))
                    )."\n";

my $BLAT_MAPPING_THRESHOLD = 80;
my $SLOP_FOR_GENOMIC_REGION = 1000;
my $SLOP_FOR_GETTING_READS = 500;
my $SMALL_INDEL_THRESHOLD = 1000;
my $REMAPPING_SCORE_THRESHOLD = 40;
my $MICROBIAL_MATCHES_MAX_FRACTION = 0.5;
my $FOOTPRINT_SIZE_MAX = 50;
my $SCORES_OUTPUT_FILE = q{};
my $BEDPEOUT = q{};
my $VIRUS_DB = q{};
my $BAC_DB_STUB = q{};
my $TMP_DIR = q{};
my $SCORE_ALG = 'ssearch36';
my $SEARCH_CORES = 1;
my $MIN_SUPPORT = 3;
GetOptions(
  'blat_mapping_threshold=i' => \$BLAT_MAPPING_THRESHOLD,
  'slop_for_genomic_region=i' => \$SLOP_FOR_GENOMIC_REGION,
  'slop_for_getting_reads=i' => \$SLOP_FOR_GETTING_READS,
  'small_indel_threshold=i' => \$SMALL_INDEL_THRESHOLD,
  'remapping_score_threshold=f' => \$REMAPPING_SCORE_THRESHOLD,
  'microbial_fraction_threshold=f' => \$MICROBIAL_MATCHES_MAX_FRACTION,
  'footprint_size_max=i' => \$FOOTPRINT_SIZE_MAX,
  'scores_output_file=s' => \$SCORES_OUTPUT_FILE,
  'virus_db=s' => \$VIRUS_DB,
  'bacterial_db_stub=s' => \$BAC_DB_STUB,
  'tmpdir=s' => \$TMP_DIR,
  'score_alg=s' => \$SCORE_ALG,
  'search_cores=i' => \$SEARCH_CORES,
  'min_support=s' => \$MIN_SUPPORT,
);

my %score_method = (emboss => \&pairwise_align_scores_emboss,
                    ssearch36 => \&pairwise_align_scores_ssearch36,
                    );
my $SCORE_SUB = $score_method{$SCORE_ALG};

my $bedpe_file = $ARGV[0];
my $bam_file = $ARGV[1];
my $fa_file = $ARGV[2];
my $bedpe_out;
$bedpe_out = $ARGV[3] if(@ARGV == 4);

if($TMP_DIR eq q{}) {
  undef $TMP_DIR;
}
else {
  make_path($TMP_DIR);
}

my $bam = Bio::DB::Sam->new(-bam => $bam_file);
my $fai = Bio::DB::Sam::Fai->load($fa_file);

my @regions;
open my $IN, '<', $bedpe_file;
while (<$IN>) {
  next if /^#/;
  chomp;
  my @F = split /\t/;
#next unless($F[6] == 3);
  push @regions, [@F];
#warn "EXIT LOOP PREMATURLY";
#last;
}
close $IN;

# For each region, do the read remapping using EMBOSS.
# At the same time aggregate reads to be mapped against the
# viral genomes.
print STDERR "Performing abnormally paired reads remapping...\n";
my (@low_end_remap_score, @high_end_remap_score, @low_end_footprint, @high_end_footprint);
my ($all_reads_fh, $all_reads_file_name) = tempfile('allreads_XXXXXX', DIR => $TMP_DIR, UNLINK=>1);
my %rg_id_of_read;

for my $i(0..$#regions) {
  my ($l, $h);
  my $r = $regions[$i];
  my $dt = DateTime->now;
  my $date = $dt->ymd .' '.$dt->hms;
  print STDERR sprintf "[%s] %d %s %s\n", $date, $i+1, (@{$r}[6,7]);
  # Deal with low end reads
  if ($r->[10] && $r->[11]) {
    $l = $r->[10];
    $h = $r->[11];
  }
  else {
    $l = ($r->[8] eq q{+} ? ($r->[1] + 60) : ($r->[2] - 100));  # These numbers are because of how Brass apparently estimates the low and high end of the ranges
    $h = ($r->[9] eq q{+} ? ($r->[4] + 60) : ($r->[5] - 100));  # These numbers are because of how Brass apparently estimates the low and high end of the ranges
  }

  if (is_small_indel($r)) {
    push @low_end_remap_score, 'small_indel';
    push @high_end_remap_score, 'small_indel';

    my @reads = collect_reads_by_region(
      $bam,
      $r->[0], $l, $r->[8],
      $r->[3], $h, $r->[9],
    );
    print_read_seqs_to_file($all_reads_fh, @reads);
    if (@reads) {
      push @low_end_footprint, get_alignment_footprint(@reads);
    }
    else {
      push @low_end_footprint, 'NA';
    }

    # Deal with high end reads
    @reads = collect_reads_by_region(
      $bam,
      $r->[3], $h, $r->[9],
      $r->[0], $l, $r->[8],
    );
    print_read_seqs_to_file($all_reads_fh, @reads);
    if (@reads) {
      push @high_end_footprint, get_alignment_footprint(@reads);
    }
    else {
      push @high_end_footprint, 'NA';
    }

    for (@reads) {
      my $read_name = $_->query->name;
      $read_name =~ s/:/_/g;
      $rg_id_of_read{$read_name} = $r->[6];
    }

    next;
  }
  else {
    # Deal with low end reads
    my $target_seq = uc(get_fai_seq($fai, $r->[3], $h-$SLOP_FOR_GENOMIC_REGION, $h+$SLOP_FOR_GENOMIC_REGION, ($r->[8] eq $r->[9] ? 1 : 0)));
    my $source_seq = uc(get_fai_seq($fai, $r->[0], $l-$SLOP_FOR_GENOMIC_REGION, $l+$SLOP_FOR_GENOMIC_REGION, 0));
    if ($target_seq =~ /N/ || $source_seq =~ /N/) {
      push @low_end_remap_score, 'ref_seq_has_N';
      push @high_end_remap_score, 'ref_seq_has_N';
      push @low_end_footprint, 'ref_seq_has_N';
      push @high_end_footprint, 'ref_seq_has_N';
      next;
    }

    my @reads = collect_reads_by_region(
      $bam,
      $r->[0], $l, $r->[8],
      $r->[3], $h, $r->[9],
    );
    print_read_seqs_to_file($all_reads_fh, @reads);
    for (@reads) {
      my $read_name = $_->query->name;
      $read_name =~ s/:/_/g;
      $rg_id_of_read{$read_name} = $r->[6];
    }
    my ($final_score, $footprint);
    if (@reads >= $MIN_SUPPORT) {
      my @low_end_score_diffs = get_remapping_score_differences(
        'low_'.$r->[6],
        \@reads,
        $target_seq,
        $source_seq,
        $r->[0], $l, $r->[8],
        $r->[3], $h, $r->[9],
      );
      $final_score = max(@low_end_score_diffs);
      $footprint = get_alignment_footprint(@reads);
    }
    else {
      $final_score = 'no_reads';
      $footprint = 'NA';
    }
    push @low_end_remap_score, $final_score;
    push @low_end_footprint, $footprint;

    # Deal with high end reads
    @reads = collect_reads_by_region(
      $bam,
      $r->[3], $h, $r->[9],
      $r->[0], $l, $r->[8],
    );
    print_read_seqs_to_file($all_reads_fh, @reads);
    for (@reads) {
      my $read_name = $_->query->name;
      $read_name =~ s/:/_/g;
      $rg_id_of_read{$read_name} = $r->[6];
    }
    $target_seq = get_fai_seq($fai, $r->[0], $l-$SLOP_FOR_GENOMIC_REGION, $l+$SLOP_FOR_GENOMIC_REGION, ($r->[8] eq $r->[9] ? 1 : 0));
    $source_seq = get_fai_seq($fai, $r->[3], $h-$SLOP_FOR_GENOMIC_REGION, $h+$SLOP_FOR_GENOMIC_REGION, 0);
    if (@reads >= $MIN_SUPPORT) {
      my @high_end_score_diffs = get_remapping_score_differences(
        'high_'.$r->[6],
        \@reads,
        $target_seq,
        $source_seq,
        $r->[3], $h, $r->[9],
        $r->[0], $l, $r->[8],
      );
      $final_score = max(@high_end_score_diffs);
      $footprint = get_alignment_footprint(@reads);
    }
    else {
      $final_score = 'no_reads';
      $footprint = 'NA';
    }
    push @high_end_remap_score, $final_score;
    push @high_end_footprint, $footprint;
  }
}

# All reads file against viral genomes
my $dt = DateTime->now;
my $date = $dt->ymd .' '.$dt->hms;
print STDERR sprintf "[%s] Running BLAT on reads against NCBI viral database...\n", $date;

my %microbe_mappings_of_rg;
my %microbe_mapping_score_of;
blat_db(\%microbe_mappings_of_rg, \%microbe_mapping_score_of, \%rg_id_of_read, $all_reads_file_name, $VIRUS_DB);


# All reads file against bacterial genomes
$date = `date`;chomp $date;
print STDERR sprintf "[%s] Running BLAT on reads against NCBI bacterial database...\n", $date;
for my $i(1..4) {
  blat_db(\%microbe_mappings_of_rg, \%microbe_mapping_score_of, \%rg_id_of_read, $all_reads_file_name, (sprintf '%s.%d.fa.2bit', $BAC_DB_STUB, $i));
}

# Print results
my $SCORES_FILE;
if ($SCORES_OUTPUT_FILE) {
  open $SCORES_FILE, '>', $SCORES_OUTPUT_FILE || die $!;
}

my $bedpe_fh;
if($bedpe_out) {
  open $bedpe_fh, '>', $bedpe_out || die $!;
}
else {
  $bedpe_fh = *STDOUT;
}

my $normalised_remap_score = $REMAPPING_SCORE_THRESHOLD * 2;
# had to double the score/penalty values as ssearch36 doesn't allow floats in gap-extension

for my $i (0..$#regions) {
  # Compile microbial mapping results
  my @microbe_mappings = values(%{$microbe_mappings_of_rg{$regions[$i]->[6]}});
  my %count_of_microbe;
  $count_of_microbe{$_}++ for @microbe_mappings;
  my @cur_microbes = sort keys(%count_of_microbe);
  my ($microbe_counts, $microbes);
  if (@cur_microbes) {
    $microbe_counts = join(',', @count_of_microbe{@cur_microbes});
    $microbes = join(',', @cur_microbes);
  }
  else {
    $microbe_counts = 0;
    $microbes = 'none';
  }

  if ($SCORES_FILE) {
    print $SCORES_FILE join(
      "\t",
      $regions[$i]->[6],
      $regions[$i]->[7],
      $low_end_remap_score[$i],
      $high_end_remap_score[$i],
      $low_end_footprint[$i],
      $high_end_footprint[$i],
      $microbe_counts,
      $microbes
    ) . "\n";
  }

  if ($low_end_remap_score[$i] eq 'no_reads') {
    print STDERR 'WARNING: no reads found from the BRM file for rearrangement ' . $regions[$i]->[6] . "!\n";
    print STDERR join("\t", @{$regions[$i]}),"\n";
  }

  # The actual filtering
  if (
    $low_end_remap_score[$i] eq 'no_reads' || $high_end_remap_score[$i] eq 'no_reads' ||
    $low_end_remap_score[$i] eq 'ref_seq_has_N' || $high_end_remap_score[$i] eq 'ref_seq_has_N' ||
    $low_end_footprint[$i] eq 'NA' || $high_end_footprint[$i] eq 'NA' ||
    ($low_end_remap_score[$i]  ne 'small_indel' && $low_end_remap_score[$i]  <= $normalised_remap_score) ||
    ($high_end_remap_score[$i] ne 'small_indel' && $high_end_remap_score[$i] <= $normalised_remap_score) ||
    sum(split /,/, $microbe_counts) >= $regions[$i]->[7] ||
    $low_end_footprint[$i] <= $FOOTPRINT_SIZE_MAX || $high_end_footprint[$i] <= $FOOTPRINT_SIZE_MAX
  ) {
    next;
  }
  print $bedpe_fh join("\t", @{$regions[$i]}) . "\n";
}

# Clean-up
close $all_reads_fh or warn $!;
#remove_tree($TMP_DIR) if(defined $TMP_DIR);
exit 0;

sub blat_db {
  my ($microbe_mappings_of_rg, $microbe_mapping_score_of, $rg_id_of_read, $reads, $db) = @_;
  my $blat = Sanger::CGP::Brass::Implement::_which('blat');
  my $blat_result = "$reads.blat";
  my $command = sprintf '%s %s %s -maxIntron=20 %s', $blat, $db, $reads, $blat_result;
  my ($stdin, $stderr, $exit) = capture { system($command); };
  die "FAILED: $command\n" if($exit != 0);
  print STDERR $stderr;
  print STDERR $stdin;
  open my $IN, '<', $blat_result || die $!;
  #$_ = <$IN> for 1..5;
  while (<$IN>) {
    next while($. < 6);
    chomp;
    my @F = split /\t/;
    next unless $F[0] >= $BLAT_MAPPING_THRESHOLD;
    if (
      !exists($microbe_mappings_of_rg->{$rg_id_of_read{$F[9]}}) ||
      !exists($microbe_mappings_of_rg->{$rg_id_of_read{$F[9]}}->{$F[9]}) ||
      $F[0] > $microbe_mapping_score_of->{$F[9]}
    ) {
      $microbe_mappings_of_rg->{$rg_id_of_read->{$F[9]}}->{$F[9]} = $F[13];
      $microbe_mapping_score_of->{$F[9]} = $F[0];
    }
  }
  close $IN;
  unlink($blat_result);
}

sub get_remapping_score_differences {

  my $identifier = shift;
  my $reads_ref = shift;
  my $target_seq = shift;
  my $source_seq = shift;

  my( $src_chr, $src_pos, $src_dir,
    $tgt_chr, $tgt_pos, $tgt_dir) = @_;

  # Source region and reads are never revcomped. Target region
  # is revcomped if strands are the same in source and target.

  my($seq_fh, $seq_fh_name) = tempfile('seq_'.$identifier.'_XXXX', DIR => $TMP_DIR, UNLINK=>1);
  print_read_seqs_to_file($seq_fh, @{$reads_ref});

  # Target region sequence
  my($target_fh, $target_fh_name) = tempfile('target_'.$identifier.'_XXXX', DIR => $TMP_DIR, UNLINK=>1);
  print($target_fh ">$tgt_chr:" . ($tgt_pos-$SLOP_FOR_GENOMIC_REGION) . '-' . ($tgt_pos+$SLOP_FOR_GENOMIC_REGION) . ':' . $tgt_dir . "\n");
  print($target_fh "$target_seq\n");
  $target_fh->flush();

  my $target_scores = $SCORE_SUB->($target_fh_name, $seq_fh_name, $TMP_DIR);

  # Source region
  my($source_fh, $source_fh_name) = tempfile('source_'.$identifier.'_XXXX', DIR => $TMP_DIR, UNLINK=>1);
  print($source_fh ">$src_chr:" . ($src_pos-$SLOP_FOR_GENOMIC_REGION) . '-' . ($src_pos+$SLOP_FOR_GENOMIC_REGION) . ':' . $src_dir . "\n");
  print($source_fh "$source_seq\n");
  $source_fh->flush();

  my $source_scores = $SCORE_SUB->($source_fh_name, $seq_fh_name, $TMP_DIR);

  # don't die here, we can live with this as the whole tree will be deleted later
  close $seq_fh or warn $!;
  close $source_fh or warn $!;
  close $target_fh or warn $!;

  my @diffs;
  for(keys %{$source_scores}) {
    $target_scores->{$_} ||= 0;
#warn "$identifier: $_: $source_scores->{$_}, $target_scores->{$_}\n";
    push @diffs, $source_scores->{$_} - $target_scores->{$_};
  }
  return @diffs;
}

sub collect_reads_by_region {
  my($bam, $chr, $pos, $dir, $mate_chr, $mate_pos, $mate_dir) = @_;
  my ($start, $end, $mate_start, $mate_end);
  # Adjust regions of expected read positions based on BEDPE intervals
  if ($dir eq '+') {
    $start = $pos - $SLOP_FOR_GETTING_READS;
    $end   = $pos + $SLOP_FOR_GETTING_READS;
  }
  else {
    $start = $pos - $SLOP_FOR_GETTING_READS;
    $end   = $pos + $SLOP_FOR_GETTING_READS;
  }
  if ($mate_dir eq '+') {
    $mate_start = $mate_pos - $SLOP_FOR_GETTING_READS;
    $mate_end   = $mate_pos + $SLOP_FOR_GETTING_READS;
  }
  else {
    $mate_start = $mate_pos - $SLOP_FOR_GETTING_READS;
    $mate_end   = $mate_pos + $SLOP_FOR_GETTING_READS;
  }

  # Change strands to Perl-Sam format
  $dir = ($dir eq '+' ? 1 : -1);
  $mate_dir = ($mate_dir eq '+' ? 1 : -1);

  my @aln = $bam->get_features_by_location($chr, $start, $end);

  @aln = grep {
    !is_supplementary($_->flag) &&
    $_->get_tag_values('FLAGS') =~ /PAIRED/ &&
    $_->get_tag_values('FLAGS') !~ /UNMAPPED/ &&
    $_->strand    eq $dir &&
    $_->mate_seq_id eq $mate_chr &&
    (
      ( $chr ne $mate_chr || $dir != $mate_dir ) ||
      ( $pos <= $mate_pos && $_->start <= $_->mate_start ) ||
      ( $pos  > $mate_pos && $_->start >= $_->mate_start )
    ) &&
    $_->mstrand   eq $mate_dir &&
    $_->mate_start  <= $mate_end &&
    $_->mate_start + $_->length - 1 >= $mate_start
  } @aln;
  return @aln;
}

sub print_read_seqs_to_file {
  my $fh = shift;
  for (@_) {
    my $read_name = $_->query->name;
    $read_name =~ s/:/_/g;
    print($fh '>', $read_name, q{ }, $_->get_tag_values('FLAGS'), ":l\n");
    print($fh $_->query->dna . "\n");
  }
  $fh->flush();
}

sub revcomp {
  my $seq = shift;
  $seq =~ tr/ACGT/TGCA/;
  return reverse $seq;
}

sub get_fai_seq {
  my $fai = shift;
  my $chr = shift;
  my $start = shift;
  my $end = shift;
  my $revcomp = shift;
  my $seq = $fai->fetch("$chr:$start-$end");
  if ($revcomp) {
    $seq = revcomp($seq);
  }
  return $seq;
}

sub pairwise_align_scores_ssearch36 {
  my $target_file = shift;
  my $seq_file = shift;
  my $tmpdir = shift;

  $target_file =~ s/^$tmpdir\///;
  $seq_file =~ s/^$tmpdir\///;

  my $ssearch36 = "cd $tmpdir;";
  $ssearch36 .= sprintf q{%s -r '2/-2' -f '-6' -g '-1' -n -m 9 -XI -3 -T %d %s %s},
                        Sanger::CGP::Brass::Implement::_which('ssearch36'),
                        $SEARCH_CORES,
                        $seq_file,
                        $target_file;

  my %scores;
  warn $ssearch36."\n";
  my $pid = open my $process, q{-|}, $ssearch36 or die 'Could not fork: '.$!;
  while(my $line = <$process>) {
    next unless($line =~ m/^>>>(.+),/);
    my @lines = ($line);
    my $name = $1;
    $line = <$process>; # blank line
    push @lines, $line;
    $line = <$process>; # '>>' line
    push @lines, $line;
    $line = <$process>; # 'sw-opt' line
    push @lines, $line;

    die qq{Failed to parse ssearch36 output, expected:\n>>>...\n\n>>...\n s-opt...\n\ngot:\n}.(join q{}, @lines)."\n" unless($lines[-2] =~ m/^>>[^>]/ && $lines[-1] =~ m/^ s-w opt:[[:space:]]+([[:digit:]]+)/);
    chomp $line;
    my ($sw_score) = $line =~ m/s-w opt:[[:space:]]+([[:digit:]]+)/;
    $scores{$name} = $sw_score if(!exists $scores{$name} || $scores{$name} < $sw_score);

  }
  close $process || die $!;
  return \%scores;
}

sub pairwise_align_scores_emboss {
  my $target_file = shift;
  my $seq_file = shift;
  my $scores_file = "$target_file.scores";

  my ($dna_fh, $dna_matrix) = tempfile('dnaMat_XXXXXX', DIR => $TMP_DIR, UNLINK=>1);
  print $dna_fh $DNA_MAT;
  $dna_fh->flush;

  my $needle = sprintf '%s -aformat score -gapo 6 -gape 1 -aseq %s -bseq %s -datafile %s -outfile %s',
                        '~yl3/programs/EMBOSS-6.6.0/bin/needle', #Sanger::CGP::Brass::Implement::_which('needle'),
                        $target_file,
                        $seq_file,
                        $dna_matrix,
                        $scores_file;

  my $exit = system($needle);
  die "Failed: $needle\n" if($exit);

  open $IN, '<', $scores_file || die $!;
  my @lines = <$IN>;
  close $IN;

  unlink($scores_file);
  close $dna_fh or warn $!;

  my %scores;
  for my $line(@lines) {
    chomp $line;
    my ($name, $score) = $line =~ m/^[+-] ([^[:space:]]+) [[:digit:]]+ \((.+)\)$/;
    next unless($name);
    $scores{$name} = $score;
  }
  return \%scores;
}

sub is_supplementary {
  return(($_[0] % 4096) >= 2048);
}

sub is_small_indel {
  my($chr1, $start1, $end1, $chr2, $start2, $end2, $dir1, $dir2) = @{$_[0]}[0..5, 8, 9];
  return(
    $chr1 eq $chr2 &&
    $dir1 eq '+' && $dir2 eq '-' && # $dir1 ne $dir2 &&
    $start2 - $end1 <= $SMALL_INDEL_THRESHOLD
  );
}

sub get_alignment_footprint {
  my @reads = @_;
  my $start = 1e9;
  my $end = 0;
  for my $a (@reads) {
    $start = $a->start if $a->start < $start;
    $end   = $a->end if $a->end > $end;
  }

  return($end - $start + 1);
}
