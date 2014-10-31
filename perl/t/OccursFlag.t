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


# for testing OccursFlag class

use File::Copy qw(copy move);
use Text::Diff;

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::BrassFilter::OccursFlag;

use Test::More 'no_plan';

use FindBin qw($Bin);

# existing entry
my $infile = $Bin.'/../testData/' . 'OccursFlag_test.in';
my $outfile = $Bin.'/../testData/' . 'OccursFlag_test.out';
my $testfile = 'OccursFlag_test.in';
my $Lfield = 23;
my $Hfield = 24;
my $bin_size = 100000;
my $within = 500;

# get the infile ready
copy $infile, $testfile;

# make a new RG object
my $OccursFlag = new Sanger::CGP::BrassFilter::OccursFlag(-infile   => $testfile,
					     -Lfield   => $Lfield,
					     -Hfield   => $Hfield,
					     -bin_size => $bin_size,
					     -within   => $within);
# process file
$OccursFlag->process();

# check the size of the outfile is the same or greater
my $infile_size = -s $testfile;
my $outfile_size = -s $outfile;

my $diff = diff "$testfile", "$outfile";

ok(defined $OccursFlag, 'object defined');

is (($OccursFlag->infile()), $testfile , "get infile");
is (($OccursFlag->Lfield()), $Lfield , "get Lfield");
is (($OccursFlag->Hfield()), $Hfield , "Hfield");
is (($OccursFlag->bin_size()), $bin_size , "bin_size");
is (($OccursFlag->within()), $within , "within");
is ($infile_size, $outfile_size , "Outfile size check");
ok(defined $diff, 'file has been modified');

if ($infile_size  == $outfile_size) { unlink $testfile; }
