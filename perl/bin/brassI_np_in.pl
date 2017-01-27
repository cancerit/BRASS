#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
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


use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => 'all';

use Capture::Tiny qw { capture };
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use File::Which qw(which);
use File::Spec;
use PCAP::Bam;
use File::Path qw(remove_tree make_path);
use FindBin qw($Bin);

const my $BAMCOLLATE2 => q{%s outputformat=sam exclude=PROPER_PAIR,UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY mapqthres=6 classes=F,F2 T=%s/bamcollate2 filename=%s};
const my $BRASS_PREP => q{%s %s -b %s.bas -np};
const my $BAMSORT => q{%s inputformat=sam verbose=0 index=1 md5=1 tmpfile=%s/bamsort md5filename=%s.md5 indexfilename=%s.bai O=%s};

die 'USAGE: <OUT_DIR> <FILE_INDEX> <FILE_1> <FILE 2>...' unless(scalar @ARGV > 2);

my $outpath = shift;
my $index = shift;
my @files = sort @ARGV;
my $bam = $files[$index-1];

my ($sample, undef) = PCAP::Bam::sample_name($bam);

my $result_path = File::Spec->catdir($outpath, $sample);
my $tmp = File::Spec->catdir($result_path, 'tmpMap');
make_path($tmp); # makes all of above
my $new_bam = File::Spec->catfile($result_path, "$sample.brm.bam");

my $command = q{};
$command .= sprintf $BAMCOLLATE2, _which('bamcollate2'), $tmp, $bam;
$command .= ' | ';
$command .= sprintf $BRASS_PREP, $^X, _which('brassI_prep_bam.pl'), $bam;
$command .= ' | ';
$command .= sprintf $BAMSORT, _which('bamsort'), $tmp, $new_bam, $new_bam, $new_bam;

warn "Starting: $command\n";

exec($command);
# replace perl with pipeline


sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}
