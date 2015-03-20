#!/usr/bin/perl

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


# filters brassI output

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::BrassFilter::BrassMarkedGroups;
use Sanger::CGP::BrassFilter::TransFlag;
use Sanger::CGP::BrassFilter::OccursFlag;
use Sanger::CGP::BrassFilter::CnFlag;
use Sanger::CGP::BrassFilter::BlatFlag;

use File::Which;
use Getopt::Long;
use Const::Fast qw(const);

my $infile = '';
my $outfile = 'brassI_filtered.out';
my $help = '';
my $debug = 0;

my $trans_only = 0;
my $occur_only = 0;
my $copyn_only = 0;
my $rblat_only = 0;

my $do_trans = 1;
my $do_occurrences = 1;
my $do_copynumber = 1;
my $do_range_blat = 1;

my $bal_trans_field = 21; # which number field of the bedpe output file the balanced translocation flag should go into
my $inv_field = 22; # which number field of the bedpe output file the inversion flag should go into
my $occL_field = 23; # which number field of the output file the L occurrences flag should go into
my $occH_field = 24; # which number field of the output file the H occurrences flag should go into
my $cn_field = 25; # which number field of the output file the near-copynumber-change flag should go into
my $blat_field = 26; # which number field of the output file the L v H blat score should go into

my $seq_depth = 30;
my $seq_depth_threshold = 25;
my $distance_threshold = 100;
my $min_tumour_count_low = 2;
my $min_tumour_count_high = 4;
my $max_normal_count = 0;
my $max_np_sample_count = 0;
my $max_np_count = 0;
my $discard_if_repeats = 0;
my $tumour = '';

# for occurence counts
my $occurs_within = 500;

# for translocations
my $bal_distance = 100000; # how far away the breakpoint coordinates can be for it to qualify as a balanced translocation
my $inv_distance = 1000; # how far away the breakpoint coordinates can be for it to qualify as an inversion

# for copynumber flagging
my $cn_within = 100000;
my $infile_ascat = '';
my $infile_ngs = '';
my $infile_bb = '';
my $cn_one_end_hit = 0;

# for blat flagging
my $ref = ''; # eg '/nfs/cancer_ref01/human/37/genome.fa';
my $blat_script = ''; # eg '/software/pubseq/bin/blat';
my $minIdentity = 95;

GetOptions( 'infile:s'               => \$infile,
            'outfile:s'              => \$outfile,
            'tumour:s'                => \$tumour,
	    'trans_only'             => \$trans_only,
	    'occurs_only'            => \$occur_only,
	    'cn_only'                => \$copyn_only,
	    'blat_only'              => \$rblat_only,
	    'seq_depth:s'            => \$seq_depth,
	    'seq_depth_threshold:s'  => \$seq_depth_threshold,
	    'distance_threshold:s'   => \$distance_threshold,
	    'min_tumour_count_low:s'  => \$min_tumour_count_low,
	    'min_tumour_count_high:s' => \$min_tumour_count_high,
	    'max_normal_count:s'     => \$max_normal_count,
	    'max_np_sample_count:s'  => \$max_np_sample_count,
	    'max_np_count:s'         => \$max_np_count,
	    'discard_if_repeats:s'   => \$discard_if_repeats,
            'bal_trans_field:s'      => \$bal_trans_field,
            'inv_field:s'            => \$inv_field,
            'occL_field:s'           => \$occL_field,
            'occH_field:s'           => \$occH_field,
            'cn_field:s'             => \$cn_field,
            'blat_field:s'           => \$blat_field,
	    'bal_distance:s'         => \$bal_distance,
	    'inv_distance:s'         => \$inv_distance,
	    'cn_distance:s'          => \$cn_within,
	    'occurs_distance:s'      => \$occurs_within,
	    'cn_one_end_hit'         => \$cn_one_end_hit,
	    'ascat:s'                => \$infile_ascat,
	    'ngs:s'                  => \$infile_ngs,
	    'bb:s'                   => \$infile_bb,
	    'ref:s'                  => \$ref,
	    'blat:s'                 => \$blat_script,
	    'minIdentity:s'          => \$minIdentity,
	    'help'                   => \$help);

# check inputs
if ($help) { usage(); }

unless ($infile) { warn "-infile not supplied\n"; usage(1); }
unless (-s $infile) { warn "-infile is empty\n"; usage(1); }
unless ($tumour)  { warn "-tumour (name of tumour) not supplied\n"; usage(1); }

# set flag constants
const my $SEQ_DEPTH => $seq_depth;
const my $SEQ_DEPTH_THRESHOLD => $seq_depth_threshold;
const my $DISTANCE_THRESHOLD => $distance_threshold;
const my $MIN_TUMOUR_COUNT_LOW => $min_tumour_count_low;
const my $MIN_TUMOUR_COUNT_HIGH => $min_tumour_count_high;
const my $MAX_NORMAL_COUNT => $max_normal_count;
const my $MAX_NP_SAMPLE_COUNT => $max_np_sample_count;
const my $MAX_NP_COUNT => $max_np_count;
const my $DISCARD_IF_REPEATS => $discard_if_repeats;

my $do_process = 1;

if ($trans_only) { $do_process = 0; $do_trans = 1; $do_occurrences = 0; $do_copynumber = 0; $do_range_blat = 0; }
if ($occur_only) { $do_process = 0; $do_trans = 0; $do_occurrences = 1; $do_copynumber = 0; $do_range_blat = 0; }
if ($copyn_only) { $do_process = 0; $do_trans = 0; $do_occurrences = 0; $do_copynumber = 1; $do_range_blat = 0; }
if ($rblat_only) { $do_process = 0; $do_trans = 0; $do_occurrences = 0; $do_copynumber = 0; $do_range_blat = 1; }

# process core information and print to outfile
if ($do_process) { process($infile, $outfile, $tumour); }

# check for translocations
if ($do_trans) { update_translocations($outfile, $bal_trans_field, $inv_field, $bal_distance, $inv_distance); }

# count the occurrences fields
if ($do_occurrences) { update_occurrences($outfile, $occL_field, $occH_field, $occurs_within); }

# check for copynumber changepoints
if ($do_copynumber)  { update_cn($outfile, $cn_field, $cn_within, $infile_ascat, $infile_ngs, $infile_bb, $cn_one_end_hit); }

# check for L v H range blat scores
if ($do_range_blat)  { update_blat($outfile, $blat_field, $blat_script, $ref, $minIdentity); }

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
sub process {
    my ($infile, $outfile, $tumour) = @_;

    unless ($outfile) { $outfile = $infile . '.bedpe'; }

    # make a summary file object
    my $file_o = new Sanger::CGP::BrassFilter::BrassMarkedGroups(-infile               => $infile,
						    -outfile              => $outfile,
						    -tumour                => $tumour,
						    -seq_depth            => $SEQ_DEPTH,
						    -seq_depth_threshold  => $SEQ_DEPTH_THRESHOLD,
						    -distance_threshold   => $DISTANCE_THRESHOLD,
						    -min_tumour_count_low  => $MIN_TUMOUR_COUNT_LOW,
						    -min_tumour_count_high => $MIN_TUMOUR_COUNT_HIGH,
						    -max_normal_count     => $MAX_NORMAL_COUNT,
						    -max_np_sample_count  => $MAX_NP_SAMPLE_COUNT,
						    -max_np_count         => $MAX_NP_COUNT,
						    -discard_if_repeats   => $DISCARD_IF_REPEATS);

    # process it
    $file_o->process();
    if ($debug) {
	my $date = `date`;
	print "finished printing to $outfile at $date";
    }
}
#------------------------------------------------------------------------------------------------#

sub update_translocations {
    my ($file, $bal_trans_field, $inv_field, $bal_distance, $inv_distance) = @_;

    unless ($file =~ /\.bedpe$/) { $file .= '.bedpe'; }

    my $trans_o = new Sanger::CGP::BrassFilter::TransFlag(-infile          => $file,
                                             -bal_distance    => $bal_distance,
                                             -inv_distance    => $inv_distance,
                                             -bal_trans_field => $bal_trans_field,
                                             -inv_field       => $inv_field);
    $trans_o->process();

    if ($debug) {
	my $date = `date`;
	print "finished trans flagging of $file at $date";
    }
}
#-----------------------------------------------------------------------------#
sub update_occurrences {
    my ($file, $occL_field, $occH_field, $within) = @_;

    unless ($file =~ /\.bedpe$/) { $file .= '.bedpe'; }

    my $occurs_o = new Sanger::CGP::BrassFilter::OccursFlag(-infile => $file,
					       -Lfield => $occL_field,
					       -Hfield => $occH_field,
                                               -within => $within );
    $occurs_o->process();

    if ($debug) {
	my $date = `date`;
	print "finished occurrences flagging of $file at $date";
    }
}
#-----------------------------------------------------------------------------#
sub update_cn {
    my ($file, $cn_field, $within, $infile_ascat, $infile_ngs, $infile_bb, $one_end_hit) = @_;

    unless ($infile_ascat || $infile_ngs || $infile_bb) {
	print "WARN: can not do copynumber changepoint flagging. No copynumber segment bed files supplied (ascat/ngs/bb)\n";
	return;
    }

    unless ($file =~ /\.bedpe$/) { $file .= '.bedpe'; }

    my $CnFlag = new Sanger::CGP::BrassFilter::CnFlag(-infile => $file,
					 -ascat  => $infile_ascat,
					 -ngs    => $infile_ngs,
					 -bb     => $infile_bb,
					 -field  => $cn_field,
					 -one_end_hit => $one_end_hit,
					 -within => $within );
    # process file
    $CnFlag->process();

    if ($debug) {
	my $date = `date`;
	print "finished CN flagging of $file at $date";
    }
}
#-----------------------------------------------------------------------------#
sub update_blat {
    my ($file, $blat_field, $blat_script, $ref, $minIdentity) = @_;

    unless ($ref && (-e "$ref")) {
	print "WARN: can not do LvH blat flagging. No valid reference file supplied\n";
	return;
    }

    if ($blat_script) {
	unless ((-e "$blat_script") && ($blat_script =~ /blat/)) {
	    my $path_to_blat = which($blat_script);
	    unless ($path_to_blat && (-e "$path_to_blat") && ($path_to_blat =~ /blat/)) {
		print "WARN: can not do LvH blat flagging. No valid blat executable supplied\n";
		return;
	    }
	}
    }
    else {
	print "WARN: can not do LvH blat flagging. No blat executable supplied\n";
	return;
    }

    unless ($file =~ /\.bedpe$/) { $file .= '.bedpe'; }

    my $BlatFlag = new Sanger::CGP::BrassFilter::BlatFlag(-infile      => $file,
					     -ref         => $ref,
					     -blat        => $blat_script,
					     -minIdentity => $minIdentity,
					     -field       => $blat_field );
    # process file
    $BlatFlag->process();

    if ($debug) {
	my $date = `date`;
	print "finished blat flagging of $file at $date";
    }
}
#-----------------------------------------------------------------------------#


sub usage {
  my $exit_code = shift;
  $exit_code ||= 0;

    print <<HERE;

brassI_filter.pl

Description - filters brassI marked groups file and outputs to another file.

options...
    -infile                : Name of the input brassI marked groups file
    -outfile               : Name of the output file (bedpe filename extension will be appended if not supplied)
    -tumour                : Name of the tumour sample

    -trans_only            : run/rerun translocation flagging of the bedpe file
    -occurrences_only      : run/rerun occurrences flagging of the bedpe file
    -cn_only               : run/rerun near-copynumber-change flagging of the bedpe file
    -blat_only             : run/rerun blat flagging of the bedpe file

    -seq_depth             : filter flag. Sequence depth for this sample. (default = 30)
    -seq_depth_threshold   : filter flag. Use min_tumour_count_high over this value and min_tumour_count_low otherwise. (default = 25)
    -distance_threshold    : filter flag. Discard rearrangements, where chrH=chrL, that do not exceed this value (default = 100)
    -min_tumour_count_low   : filter flag. Discard rearrangements which do not reach this number of reads in any of the tumour or metastatic samples involved
                             (low seq_depth) (default = 2)
    -min_tumour_count_high  : filter flag. Discard rearrangements which do not reach this number of reads in any of the tumour or metastatic samples involved
                             (above seq_depth threshold) (default = 4)
    -max_normal_count      : filter flag. Discard rearrangements which have more than this number of reads in the matched normal (default = 0)
    -max_np_sample_count   : filter flag. Discard rearrangements which have more than this number of unmatched normal panel samples with reads  (default = 0)
    -max_np_count          : filter flag. Discard rearrangements which have more than this number of reads in the unmatched normal panel samples (default = 0)
    -discard_if_repeats    : filter flag. Discard rearrangements which are associated with known repeats (default = 0)

    -bal_trans_field       : which number field of the bedpe output file the balanced translocation flag should go into (default = 21)
    -inv_field             : which number field of the bedpe output file the inversion flag should go into (default = 22)
    -occL_field            : which number field of the bedpe output file the L occurrences flag should go into (default = 23)
    -occH_field            : which number field of the bedpe output file the H occurrences flag should go into (default = 24)
    -cn_field              : which number field of the bedpe output file the near-copynumber-change flag should go into (default = 25)
    -blat_field            : which number field of the bedpe output file the L v H blat score should go into (default = 26)

    -occurs_distance       : how far away breakpoint ends in 2 different rearrangements can be, to be declared a similar coordinate (default = 500)

    -bal_distance          : how far away the breakpoint coordinates for 2 different rearrangements can be, to be declared a balanced translocation (default = 100000)
    -inv_distance          : how far away the breakpoint coordinates for 2 different rearrangements can be, to be declared a inversion (default = 1000)

    -cn_distance           : how far away a copynumber chagepoint can be from a rearrangement for it to qualify as a changepoint hit (default = 100000)
    -cn_one_end_hit        : If -cn_one_end_hit is present, call a hit even if only one end of the rearrangement is near a changepoint (default = not set)
    -ascat                 : copynumber (cn) segments summary file for this sample - ASCAT. (optional)
                             Line Format: unused,chr,start,end,normal_total_cn(optional),normal_minor_cn(optional),tumour_total_cn,tumour_minor_cn
    -ngs                   : copynumber (cn) segments summary file for this sample - NGS. (optional)
                             Line Format: unused,chr,start,end,normal_total_cn(optional),normal_minor_cn(optional),tumour_total_cn,tumour_minor_cn
    -bb                    : copynumber (cn) segments summary file for this sample - Battenberg. (optional)
                             Line Format: unused,chr,start,end,normal_total_cn(optional),normal_minor_cn(optional),tumour_total_cn,tumour_minor_cn

    -ref                   : Blat of breakpoint range L against range H - fasta format Reference file (fai index file also present) to retrieve breakpoint range sequence from.
    -blat                  : Blat of breakpoint range L against range H - blat script to use (default = blat)
    -minIdentity           : Blat of breakpoint range L against range H - minimum identity value to supply to blat (default = 95)

    -help                  : Print this message


example...

brassI_filter.pl -infile file.rg.marked -outfile test_brass.out -tumour AB1234 -ascat AB1234_ascat.summary.csv  -ngs AB1234_ngs.summary.csv -ref /nfs/cancer_ref01/human/37/genome.fa -blat "/software/pubseq/bin/blat" -cn_one_end_hit


Default output column order (tab delimited) and description:

chr1             - lower breakpoint chromosome
start1           - lower breakpoint range start (zero referenced)
end1             - lower breakpoint range end (zero referenced)
chr2             - higher breakpoint chromosome
start2           - higher breakpoint range start
end2             - higher breakpoint range end
id/name          - name or id of this rearrangement
brass_score      - brassI score (total number of reads)
strand1          - lower breakpoint strand
strand2          - higher breakpoint strand
repeats          - known repeats detected by brassI
np_sample_count  - how may normal panel samples have this rearrangement
tumour_count      - how many reads represent this rearrangement in the tumour
normal_count     - how many reads represent this rearrangement in the normal
np_count         - how many reads represent this rearrangement in the normal panel samples
bkpt_distance    - how far apart the breakpoints are (if on the same chromosome)
sample           - name of sample
sample_type      - type of sample T(umour), N(ormal), M(etastasis), N(ormal)P(anel)
names            - list of readnames that came from the stated sample
count            - number of reads
bal_trans        - list of ids/names for other rearrangement entries that may pair with this rearrangement to form a balanced translocation
inv              - list of ids/names for other rearrangement entries that may pair with this rearrangement to form an inversion
occL             - how many time the lower (1) breakpoint appears in this dataset
occH             - how many time the higher (2) breakpoint appears in this dataset
copynumber_flag  - Summary of copynumber changepoints that lie close to the Lower and Higher breakpoints of this rearrangement.
                   Format: Lower_Higher. A = Ascat, N = Picnic NGS, B = Battenberg. eg A_ANB.
range_blat       - blat score from aligning lower_breakpoint_sequence to higher_breakpoint_sequence




Author : las


HERE

exit($exit_code);
}
