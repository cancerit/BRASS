#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL=>'all';
use autodie;
use Capture::Tiny qw(capture);
use Scalar::Util qw(looks_like_number);
use Const::Fast qw(const);
use File::Which qw(which);
use FindBin qw($Bin);

use PCAP::Bam;

const my $PROPER_INTERSECT => q{bash -c 'set -o pipefail; samtools view -ub -f 2 -F 3852 -q 1 %s %s | bedtools intersect -loj -a %s -b stdin -sorted -c > %s'};
const my $FOLDBK_INTERSECT => q{bash -c 'set -o pipefail; (samtools view -H %s; samtools view -F 3854 %s %s | %s %s) | samtools view -Sbu - | bedtools intersect -loj -a %s -b stdin -sorted -c > %s'};

if(scalar @ARGV < 3) {
  die "USAGE: ./compute_coverage.pl gc_windows.bed[.gz] in.bam out_path [chr_name]\n";
}

my ($windows, $bam, $out) = @ARGV;
my $chr_name;
if(scalar @ARGV == 4) {
  $chr_name = $ARGV[3];
}

# input checks
die "ERROR: gc_windows.bed[gz] file must exist with non-zero size\n" unless(-e $windows && -s _);
die "ERROR: in.bam file must exist with non-zero size\n" unless(-e $bam && -s _);

my $sample_name = sanitised_sample_from_bam($bam);

my $cn_bed_file = qq{$out/${sample_name}.ngscn.bed.gz};
my $cn_fb_file = qq{$out/${sample_name}.ngscn.fb_reads.bed.gz};
my $new_windows;
if(defined $chr_name) {
  # need to get windows for just this chr and write to a tmp file
  # now need to capture all windows from original window file by the chr_name
  $new_windows = "$out/windows_$chr_name.bed";
  run_ext(qq{zgrep -wE '^$chr_name' $windows > $new_windows});
  $cn_bed_file = qq{$out/${sample_name}.$chr_name.ngscn.bed};
  $cn_fb_file = qq{$out/${sample_name}.$chr_name.ngscn.fb_reads.bed};
}
run_ext(sprintf $PROPER_INTERSECT, $bam, $chr_name, ($new_windows || $windows), $cn_bed_file);
run_ext(sprintf $FOLDBK_INTERSECT, $bam, $bam, $chr_name, $^X, _which('brass_foldback_reads.pl'), $cn_bed_file, $cn_fb_file);
unlink $new_windows if(-e $new_windows);

sub run_ext {
  my ($command) = @_;
  warn "Running: $command\n";
  my ($stdout, $stderr, $exit) = capture{ system($command); };
  die "ERROR: Failed to execute $command\n\tERROR: $stderr\n\tEXIT CODE: $exit\n" if($exit);
  warn "Exit Code: $exit\n";
  return $stdout;
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_\-.]/_/ig; # sanitise sample name
  return $sample;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
  return $path;
}
