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

my ($full_bed_pe, $rg_patt_pe, $out_file) = @ARGV;

my %ids;
open my $FILT, '<', $rg_patt_pe || die $!;
while(my $line = <$FILT>) {
  my ($id, $low, $high) = (split /\t/, $line)[6,12,13];
  $ids{$id} = [$low, $high];
}
close $FILT;

my $ofh = *STDOUT;
if(defined $out_file) {
  open $ofh, '>', $out_file || die $!;
}

open my $MAIN, '<', $full_bed_pe || die $!;
while(my $line = <$MAIN>) {
  if($line =~ m/^#/) {
    print $line;
    next;
  }

  my @F = split /\t/, $line;
  my $id = $F[6];
  next unless(exists $ids{$id});
  $F[1] = $ids{$id}->[0]-1;
  $F[2] = $ids{$id}->[0];
  $F[4] = $ids{$id}->[1]-1;
  $F[5] = $ids{$id}->[1];
  print join("\t", @F);
}
close $MAIN;

close $ofh if(defined $out_file);
