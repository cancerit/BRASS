package Bio::Brass::Alignment;

########## LICENCE ##########
# Copyright (c) 2014,2015 Genome Research Ltd.
#
# Author: Cancer Genome Project <cgpit@sanger.ac.uk>
#
# This file is part of BRASS.
#
# BRASS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’."
########## LICENCE ##########

use strict;
use warnings;

use base qw(Exporter);
our @EXPORT = qw(map_sequences);

use IPC::Run qw(run start new_chunker timeout);
use File::Temp qw(tempfile);
use Bio::Brass::Location;
use Bio::Brass qw($VERSION);

=head1 NAME

Bio::Brass::Alignment - represent alignments as Sugar strings

=head1 SYNOPSIS

    use Bio::Brass::Alignment;

    my @mappings = map_sequences($options, @sequences);

    $aln->trim();

=head1 DESCRIPTION

=cut

use constant { _Q => 0, _T => 1, _SCORE => 2 };
use constant { _POS => 0, _LIM => 1, _SEQ => 2 }; # FIXME import from Location

use overload
    '""' => \&to_string,
    'cmp' => \&compare_by_score,
    '<=>' => \&compare_by_query_location;

=head2 new

    $aln = Bio::Brass::Alignment->new($q_locn, $t_locn, $score);

=cut

sub new {
    my $class = shift;
    return bless [@_], $class;
}

sub to_string {
    my ($self) = @_;
    return "{q:$$self[_Q], t:$$self[_T], score $$self[_SCORE]}";
}

=head2 Comparisons

    $signum = $aln1 cmp $aln2;  # Equivalent to using by_score
    @alns = sort Bio::Brass::Alignment::by_score @alns;
    @alns = sort Bio::Brass::Alignment::by_query_location @alns;

Spaceship operators returning -1, 0, or +1 according to the ordering between
the two alignments.  With C<cmp> or C<compare_by_score()>, alignments are
ordered by score (highest first) and then by target genomic location.
With C<E<lt>=E<gt>> or C<compare_by_query_location()>, they are ordered only by query location.

=cut

sub by_score {
    my $diff = $$b[_SCORE] <=> $$a[_SCORE];
    $diff = $$a[_T]->name cmp $$b[_T]->name if $diff == 0;
    $diff = $$a[_T] cmp $$b[_T] if $diff == 0;
    return $diff;
}

sub compare_by_score {
    my ($a, $b, $reversed) = @_;
    my $diff = $$b[_SCORE] <=> $$a[_SCORE];
    $diff = $$a[_T]->name  cmp $$b[_T]->name  if $diff == 0;
    $diff = $$a[_T] cmp $$b[_T] if $diff == 0;
    return $reversed? -$diff : +$diff;
}

sub by_query_location {
#print STDERR "Args: @_\n";
#print STDERR "1:$a: 2:$$a[_Q]:\n";
#print STDERR "1:$b: 2:$$b[_Q]:\n";
    return $$a[_Q] cmp $$b[_Q];
}

sub compare_by_query_location {
    my ($a, $b, $reversed) = @_;
#print STDERR "compare_by_query_location($a, $b)\n";
    my $diff = $$a[_Q] cmp $$b[_Q];
    return $reversed? -$diff : +$diff;
}

sub query  { return ${$_[0]}[_Q] }
sub target { return ${$_[0]}[_T] }
sub score  { return ${$_[0]}[_SCORE] }

sub trim {
    my ($self, $left, $right) = @_;
    $$self[_Q]->trim($left, $right);
    $$self[_T]->trim($left, $right);
}

# jwh basically the point at which the two sequences differ continuosly for more than the given threshold
# return the current position.
sub _trim_end {
    my ($query, $qbase, $target, $tbase, $dir, $len, $threshold) = @_;
    my ($startpos, $overhang);

    my $i = 0;
    do {
      $startpos = $i;

      $i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
			               eq substr($target, $tbase + $dir * $i, 1);

      my $matched_count = $i - $startpos;
      my $mismatchpos = $i;

      $i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
                           ne substr($target, $tbase + $dir * $i, 1);

      my $mismatched_count = $i - $mismatchpos;
      $overhang = $matched_count - $mismatched_count;

      #print STDERR "_trim_end($dir): $overhang\n";
    } until $i == $len || $overhang >= $threshold;

    return $startpos;
}

sub trim_mismatches {
    my ($self, $ends, $threshold) = @_;
    $ends = 'LR' unless defined $ends;
    $threshold = 6 unless defined $threshold;

    if ($ends =~ /[Ll]/) {
	my $trim = _trim_end($$self[_Q][_SEQ], $$self[_Q][_POS],
			     $$self[_T][_SEQ], $$self[_T][_POS],
			     +1, $$self[_T]->length, $threshold);
	$$self[_Q][_POS] += $trim;
	$$self[_T][_POS] += $trim;
    }

    if ($ends =~ /[Rr]/) {
	my $trim = _trim_end($$self[_Q][_SEQ], $$self[_Q][_LIM] - 1,
			     $$self[_T][_SEQ], $$self[_T][_LIM] - 1,
			     -1, $$self[_T]->length, $threshold);
	$$self[_Q][_LIM] -= $trim;
	$$self[_T][_LIM] -= $trim;
    }

    return $self;
}

=head2 untrim

    $aln->untrim();

Undoes any trimming and other cookery, returning $aln to its original
boundaries.

=cut

sub untrim {
    my ($self) = @_;
    $$self[_Q]->untrim();
    $$self[_T]->untrim();
}

sub revcomp {
    my ($self) = @_;
    $$self[_Q]->revcomp();
    $$self[_T]->revcomp();
}

#use List::Util qw(min max);

# Given a pair of integers, returns the greater one.
#sub max {
#    my ($a, $b) = @_;
#    return ($a > $b)? $a : $b;
#}

=head2 map_sequences

    @mappings = map_sequences($options, @sequences);

Blah blah.

=cut

our $exonerate_target;
our $ref_fai;
our $max_iteration_attemps = 5;
our $itteration_wait = 5;

sub map_sequences {
	my($working_dir,$exonerate_options,$regions,@contigs) = @_;

    $exonerate_options = q{} unless defined $exonerate_options;
    my $sort_function = undef;

#    local $_;
#    die __PACKAGE__."::map_sequences(): \$exonerate_target not set\n"
#      unless defined $exonerate_target;

    die __PACKAGE__."::map_sequences(): \$ref_fai not set\n"
      unless defined $ref_fai;

    my @mapping = ();

    my $i = 0;

    my ($FASTA, $queryfname) = tempfile('brassAssExonerateXXXXXXX', DIR => $working_dir, SUFFIX => '.fa');
    foreach my $contig (@contigs) {
        if (ref $contig) { print $FASTA ">$i\n", $contig->seq(), "\n" or die __PACKAGE__."::map_sequences(): write $FASTA failed: $!\n";}
        else { print $FASTA ">$i\n$contig\n" or die __PACKAGE__."::map_sequences(): write $FASTA failed: $!\n";}
        push @mapping, [];

        $i++;
    }
    close $FASTA or die __PACKAGE__."::map_sequences(): close $FASTA failed: $!\n";

    my ($FASTA_ref, $reffname) = tempfile("brassAssExonerateRefXXXXXXX", DIR => $working_dir, SUFFIX => '.fa');
    foreach my $region (@$regions) {
        print $FASTA_ref ">$region\n" or die __PACKAGE__."::map_sequences(): write $FASTA_ref failed: $!\n";
    	print $FASTA_ref $ref_fai->fetch($region)."\n" or die __PACKAGE__."::map_sequences(): write $FASTA_ref failed: $!\n";
    }
    close $FASTA_ref or die __PACKAGE__."::map_sequences(): close $FASTA_ref failed: $!\n";

    #    open my $EXONERATE,
    #	"exonerate --verbose 0 --showalignment no --showvulgar no " .
    #	"--querytype dna --targettype dna $exonerate_options " .
    #	"--ryo '%qas*Q %qab %qae %qS\n%tas*T %tab %tae %tS %ti\n** %qi %s\n' ".
    #	"--query $queryfname --target $exonerate_target|"
    #	or die __PACKAGE__."::map_sequences(): can't spawn exonerate: $!\n";

    my @cmd = (
          'exonerate',
          '--verbose','0',
          '--showalignment','no',
          '--showvulgar','no',
          '--showcigar','no',
          '--querytype','dna',
          '--targettype','dna',
          '--ryo',q{%qas*Q %qab %qae %qS\n%tas*T %tab %tae %tS %ti\n** %qi %s\n},
          '--query',$queryfname,
          '--target',$reffname, #$exonerate_target,
          split(' ',$exonerate_options),
    );

    _run_exonerate(\@cmd,\@mapping,$working_dir,1);

    unlink $queryfname;
    unlink $reffname;

    if (defined $sort_function) {
        @$_ = sort $sort_function @$_ foreach @mapping;
    }

    return @mapping;
}

sub _run_exonerate{
    my($cmd,$mapping,$working_dir,$itteration) = @_;

	my $err = '';
  	my $seq = '';
    my @locations = ();
    my $mappings = 0; ## no itteration if fails while streaming data.
	my $fh;
  warn "$cmd\n" if($itteration == 1);
	eval{
        run($cmd, '>', new_chunker, sub{
            my ($line) = @_;
            chomp $line;

            if (substr($line, 0, 1) ne '*') {
                $seq .= $line;
            }
            elsif (substr($line, 0, 2) ne '**') {
                #print STDERR "EXON: {$line}\n";
                #T %tab %tae %tS %ti
                my ($type, $begin, $end, $strand, $reg) = split /\s+/, $line, 5;

                ## need to add the actual genomic coords to the target
                if($type eq '*T'){
                	my($chr,$srt,$stp) = $reg =~ /(.+):(\d+)-(\d+)/;
                	$reg = $chr if defined $chr;
                	$begin += $srt - 1 if defined $srt;
                	$end += $srt - 1 if defined $stp;
                }

                push @locations, Bio::Brass::Location->new($seq, $begin, $end, $strand, $reg);
                $seq = "";
                #my $lenfoo = $locations[-1]->length;
                #print STDERR "LOCN: $locations[-1]; length: $lenfoo\n";
            }
            else {
                my (undef, $query, $score) = split /\s+/, $line, 3;
                my $aln = Bio::Brass::Alignment->new(@locations, $score);
                $aln->revcomp() if $locations[0]->strand eq '-';

                ## The following prevents multiple alignments to the same place when overlapping regions are used to align
                ## against.
                ## There is the possibility a contig could map to both regions but only partially to one (in cases where
                ## the serach space of the other contig parially overlaps with the current contig).
                ## If we find an exact or overlapping mapping in our current loccations replace it with the one being read
                ## if its score is equal or higher.
                my $cur_start = $aln->target->start;
                my $cur_end = $aln->target->end;
                ($cur_start,$cur_end) = ($cur_end,$cur_start) if ($cur_end < $cur_start);
                my $replace = -1;
                for(my $i=0;$i<scalar(@{$mapping->[$query]});$i++){
                    my $existing_aln = $mapping->[$query]->[$i];

                    if($existing_aln->target->name eq $aln->target->name){
                        my $ex_start = $existing_aln->target->start;
                        my $ex_end = $existing_aln->target->end;
                        ($ex_start,$ex_end) = ($ex_end,$ex_start) if ($ex_end < $ex_start);

                	    if( ($ex_start >= $cur_start && $ex_start <= $cur_end)
                	      || ($ex_end <= $cur_end && $ex_end >= $cur_start)
                	      || ($cur_start < $ex_start && $cur_end > $ex_end)
                	      ){
                	        if($score >= $existing_aln->score){
                	            $replace = $i;
                	            last;
                	        }else{
                	        	#if it overlaps but not any better in score skip it!
                	        	$replace = -2;
                	        }
                	    }
                    }
                }

                if($replace > -1){
                    $mapping->[$query]->[$replace] = $aln;
                }elsif($replace == -1){
                    push @{$mapping->[$query]}, $aln;
                }

                #print STDERR "ALN:  $aln\n";
                @locations = ();
                $mappings++;
            }
            return 1;
        }, '2>', sub{
            $err .= $_[0];
        },timeout( 200 )) or die "error code: $?, msg: $err";

    };if($@){
    	close $fh if(defined $fh);
    	if($@ =~ /IPC::Run: timeout on timer/g){
    	    warn __PACKAGE__."::_run_exonerate(): SKIPPING on itteration: $itteration (of $max_iteration_attemps), $@";
    	}else{

            if($itteration >= $max_iteration_attemps || $mappings){
    	        die __PACKAGE__."::_run_exonerate(): Exonerate error: killed at itteration: $itteration (of $max_iteration_attemps), $@";
            }
            warn __PACKAGE__."::_run_exonerate(): Exonerate error: itteration: $itteration (of $max_iteration_attemps), $@";
            eval{
               sleep($itteration_wait);
            };if($@){
                die __PACKAGE__."::_run_exonerate(): Sleep error: itteration: $itteration (of $max_iteration_attemps), msg $@" unless $@ eq "Alarm!\n";
            }
            _run_exonerate($cmd,$mapping,$working_dir,++$itteration);
    	}
    }
}

=head1 AUTHOR

John Marshall E<lt>jm18@sanger.ac.ukE<gt>

=cut

1;
__END__

=head2 trim_alignment

    $subsugar = trim_alignment($seqtext, $subsugar, $refseqtext);
    $subsugar = trim_alignment($st, $ss, $rst, $ends, $threshold);

Blah blah.

=cut

sub _trim_end {
    my ($query, $qbase, $target, $tbase, $dir, $len, $threshold) = @_;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
			     eq substr($target, $tbase + $dir * $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
			     ne substr($target, $tbase + $dir * $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $threshold;

    return $startpos;
}

sub trim_alignment {
    my ($query, $subsugar, $subtarget, $ends, $threshold) = @_;
    $ends = 'LR' unless defined $ends;
    $threshold = 6 unless defined $threshold;

    my ($qpos, $qlim, $qstrand, $tname, $tpos, $tlim, $tstrand, $score) =
	split /\s+/, $subsugar;

    if ($ends =~ /[Ll]/) {
	my $trim = _trim_end($query, $qpos, $subtarget, 0, +1,
			     length $subtarget, $threshold);
	$qpos += $trim;
	$tpos += $trim;
    }

    if ($ends =~ /[Rr]/) {
	my $trim = _trim_end($query, $qlim - 1, $subtarget, -1, -1,
			     length $subtarget, $threshold);
	$qlim -= $trim;
	$tlim -= $trim;
    }

    return "$qpos $qlim $qstrand $tname $tpos $tlim $tstrand $score";
}

=head2 trim

    $aln->trim(...);

=cut

sub trim {
    my ($self, $cooker, $threshold) = @_;

    $threshold = 6 unless defined $threshold;
    $cooker = "basic" unless defined $cooker;

    local $_ = "$cooker-$threshold";
    if (exists $$self[COOKER]) { $$self[COOKER] .= " $_" }
    else { push @$self, $_, @$self[QPOS, QLIM, TPOS, TLIM] }

    $_ = _trim_end($$self[RAW_QSEQ], $$self[QPOS],
		   $$self[RAW_TSUBSEQ], $$self[TPOS] - $$self[RAW_TPOS],
		   +1, $$self[QLIM] - $$self[QPOS], $threshold);
print "TRIM $self: pos += $_\n" if $_ > 0 && $$self[TSTRAND] eq '-';
    $$self[QPOS] += $_;
    $$self[TPOS] += $_;

    $_ = _trim_end($$self[RAW_QSEQ], $$self[QLIM] - 1,
		   $$self[RAW_TSUBSEQ], -1 - ($$self[RAW_TLIM] - $$self[TLIM]),
		   -1, $$self[QLIM] - $$self[QPOS], $threshold);
print "TRIM $self: lim -= $_\n" if $_ > 0 && $$self[TSTRAND] eq '-';
    $$self[QLIM] -= $_;
    $$self[TLIM] -= $_;
}

=head2 trim_alignment

    Bio::Brass::Alignment::trim_threshold = 6;
    @subsugar = trim_alignment_left ($seqtext, $refseqtext, @subsugar);
    @subsugar = trim_alignment_right($seqtext, $refseqtext, @subsugar);
    @subsugar = trim_alignment      ($seqtext, $refseqtext, @subsugar);

Blah blah.

=cut

our $trim_threshold = 6;

sub trim_alignment_left {
    my ($query, $subtarget, @subsugar) = @_;
    my $qpos = $subsugar[QPOS];
    my $len = length $subtarget;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query, $qpos + $i, 1)
			     eq substr($subtarget, $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query, $qpos + $i, 1)
			     ne substr($subtarget, $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $trim_threshold;

    $subsugar[QPOS] += $startpos;
    $subsugar[TPOS] += $startpos;
    return @subsugar;
}

#	my $trim = _trim_end($query, $qlim - 1, $subtarget, -1, -1,
#			     length $subtarget, $threshold);
sub trim_alignment_right {
    my ($query, $subtarget, @subsugar) = @_;
    my $qend = $subsugar[QLIM] - 1;
    my $len = length $subtarget;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query, $qend - $i, 1)
			     eq substr($subtarget, -1 - $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query, $qend - $i, 1)
			     ne substr($subtarget, -1 - $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $trim_threshold;

    $subsugar[QLIM] -= $startpos;
    $subsugar[TLIM] -= $startpos;
    return @subsugar;
}

sub trim_alignment_bork {
    my ($query, $subtarget, @subsugar) = @_;

    @subsugar = trim_alignment_left ($query, $subtarget, @subsugar);
    @subsugar = trim_alignment_right($query, $subtarget, @subsugar);
    return @subsugar;

    return trim_alignment_left($query, $subtarget,
		    trim_alignment_right($query, $subtarget, @subsugar));
}

