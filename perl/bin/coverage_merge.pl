#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use autodie;
use warnings FATAL => 'all';
use Cwd;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);
use List::Util qw(first);

const my $FILE_PAIR => '%s.%s.ngscn.bed';
const my $FILE_FOLD => '%s.%s.ngscn.fb_reads.bed';

die "Usage: genome.fa.fai sample_name indir include_list" unless(scalar @ARGV >= 3);

my ($fai, $sample, $indir, $include) = @ARGV;
die "ERROR: *.fai file must exist with non-zero size\n" unless(-e $fai && -s _);
die "ERROR: indir must exist\n" unless(-e $indir && -d _);
die "ERROR: include_list must be provided (csv of valid chrs to use)\n" unless(defined $include);

my @inc_list = split /,/, $include;

my @chr_order;
open my $FAI_IN, '<', $fai;
while(<$FAI_IN>) {
  my $chr = (split /\t/)[0];
  push @chr_order, $chr if(first {$chr eq $_} @inc_list);
}
close $FAI_IN;

my $final_pair = "$indir/$sample.ngscn.bed.gz";
my $final_fold = "$indir/$sample.ngscn.fb_reads.bed.gz";

unlink $final_pair if(-e $final_pair);
unlink $final_fold if(-e $final_fold);

my $init_dir = getcwd;

my $err_code = 0;
try {
  chdir $indir;
  cat_to_gzip($FILE_PAIR, $final_pair, $sample, \@chr_order);
  cat_to_gzip($FILE_FOLD, $final_fold, $sample, \@chr_order);
} catch {
  if($_) {
    warn $_;
    $err_code = 1;
  }
} finally {
  chdir $init_dir;
};

exit $err_code;

sub cat_to_gzip {
  my ($format, $outfile, $sample, $chrs) = @_;
  my @args;
  for my $chr(@{$chrs}) {
    push @args, sprintf $format, $sample, $chr;
    die "Expected file missing $indir/$args[-1]\n" unless(-e $args[-1]);
  }
  my $command = qq{bash -c 'set -o pipefail; cat @args | gzip -c > $outfile'};
  warn $command."\n";
  system($command) == 0 or die "Failed to merge files to $outfile: $!\n";
}
