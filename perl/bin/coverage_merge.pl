#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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
use autodie;
use warnings FATAL => 'all';
use Cwd;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);
use List::Util qw(first);

const my $FILE_PAIR => '%s.%s.ngscn.bed';
const my $FILE_FOLD => '%s.%s.ngscn.fb_reads.bed';

die "Usage: genome.fa.fai sample_name indir [exclude_list]" unless(scalar @ARGV >= 3);

my ($fai, $sample, $indir, $exclude) = @ARGV;
die "ERROR: *.fai file must exist with non-zero size\n" unless(-e $fai && -s _);
die "ERROR: indir must exist\n" unless(-e $indir && -d _);

my @chr_order;
open my $FAI_IN, '<', $fai;
while(<$FAI_IN>) {
  push @chr_order, (split /\t/)[0];
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
  my $exclude_list = exclude_patterns($exclude);
  cat_to_gzip($FILE_PAIR, $final_pair, $sample, \@chr_order, $exclude_list);
  cat_to_gzip($FILE_FOLD, $final_fold, $sample, \@chr_order, $exclude_list);
} catch {
  if($_) {
    warn $_;
    $err_code = 1;
  }
} finally {
  chdir $init_dir;
};

exit $err_code;

sub exclude_patterns {
  my $patt = shift;
  my @exclude;
  return \@exclude unless($patt);
  @exclude = split /,/, $patt;
  my @exclude_patt;
  for my $ex(@exclude) {
    $ex =~ s/%/.+/;
    push @exclude_patt, $ex;
  }
  return \@exclude;
}

sub cat_to_gzip {
  my ($format, $outfile, $sample, $chrs, $exclude_list) = @_;
  my @args;
  for my $chr(@{$chrs}) {
    next if(first { $chr =~ m/^$_$/ } @{$exclude_list});
    push @args, sprintf $format, $sample, $chr;
    die "Expected file missing $indir/$args[-1]\n" unless(-e $args[-1]);
  }
  my $command = qq{bash -c 'set -o pipefail; cat @args | gzip -c > $outfile'};
  warn $command."\n";
  system($command) == 0 or die "Failed to merge files to $outfile: $!\n";
}
