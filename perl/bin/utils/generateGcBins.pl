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
use Const::Fast qw(const);
use Bio::DB::HTS::Faidx;

const my $BIN_SIZE => 500;

die "USAGE: generateGcBins.pl genome.fa > gcBins.bed\n" unless(@ARGV);

my $ref = shift @ARGV;
my $index = Bio::DB::HTS::Faidx->new($ref);

my @seq_ids = $index->get_all_sequence_ids();
for my $seq_id(@seq_ids) {
  my $len = $index->length($seq_id);
  warn "Processing seq: $seq_id ($len b.p.)\n";
  my $low = 0;
  while($low < $len) {
    my $high = $low+$BIN_SIZE;
    $high = $len if($high > $len);
    my $seq = $index->get_sequence_no_length(sprintf "%s:%d-%d", $seq_id, $low+1, $high);
    $seq =~ s/N//g;
    my $s_len = length $seq;
    my $gc = $seq =~ tr/GC//;
    if($gc == 0) {
      $gc = 'NA'
    }
    else {
      $gc = sprintf '%.6f', $gc/$s_len;
    }
    printf "%s\t%d\t%d\t%d\t%s\n", $seq_id, $low, $high, $s_len, $gc;
    $low+= $BIN_SIZE;
  }
}
