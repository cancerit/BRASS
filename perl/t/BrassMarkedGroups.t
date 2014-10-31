#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd.
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
########## LICENCE ##########


# for testing BrassMarkedGroups class

use File::Copy qw(copy move);
use Text::Diff;

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::BrassFilter::BrassMarkedGroups;

use Test::More 'no_plan';

use FindBin qw($Bin);

# existing entry
my $infile = $Bin.'/../testData/' . 'BrassMarkedGroups_test.in';
my $outfile = $Bin.'/../testData/' . 'BrassMarkedGroups_test.out.bedpe';
my $testfile = 'BrassMarkedGroups_test.in';
my $testoutfile = 'BrassMarkedGroups_test.out.bedpe';

my $tumour = 'AB1234';
my $min_tumour_count_high = 3;
my $min_tumour_count_low = 1;
my $max_normal_count = 3;
my $max_np_count = 5;
my $max_np_sample_count = 2;
my $distance_threshold = 1000;
my $seq_depth_threshold = 30;
my $seq_depth = 20;
my $discard_if_repeats = 0;


# get the infile ready
copy $infile, $testfile;

# make a new object
my $BrassMarkedGroups = new Sanger::CGP::BrassFilter::BrassMarkedGroups(-infile              => $testfile,
							   -outfile             => $testoutfile,
							   -tumour               => $tumour,
							   -min_tumour_count_high=> $min_tumour_count_high,
							   -min_tumour_count_low => $min_tumour_count_low,
							   -max_normal_count    => $max_normal_count,
							   -max_np_count        => $max_np_count,
							   -max_np_sample_count => $max_np_sample_count,
							   -distance_threshold  => $distance_threshold,
							   -seq_depth_threshold => $seq_depth_threshold,
							   -seq_depth           => $seq_depth,
							   -discard_if_repeats  => $discard_if_repeats);
# process file
$BrassMarkedGroups->process();

my $diff = diff "$testoutfile", "$outfile";

ok(defined $BrassMarkedGroups, 'object defined');

is (($BrassMarkedGroups->infile()), $testfile , "get infile");
is (($BrassMarkedGroups->tumour()), $tumour , "get tumour");
is (($BrassMarkedGroups->min_tumour_count_high()), $min_tumour_count_high , "get min_tumour_count_high");
is (($BrassMarkedGroups->min_tumour_count_low()), $min_tumour_count_low , "get min_tumour_count_low");
is (($BrassMarkedGroups->max_normal_count()), $max_normal_count , "get max_normal_count");
is (($BrassMarkedGroups->max_np_count()), $max_np_count , "get max_np_count");
is (($BrassMarkedGroups->max_np_sample_count ()), $max_np_sample_count  , "get max_np_sample_count ");
is (($BrassMarkedGroups->distance_threshold()), $distance_threshold , "get distance_threshold");
is (($BrassMarkedGroups->seq_depth_threshold()), $seq_depth_threshold , "get seq_depth_threshold");
is (($BrassMarkedGroups->seq_depth()), $seq_depth , "get seq_depth");
is (($BrassMarkedGroups->discard_if_repeats()), $discard_if_repeats , "get discard_if_repeats");
ok(!$diff, 'correct file created');

unless ($diff) { unlink $testfile; unlink $testoutfile; }
