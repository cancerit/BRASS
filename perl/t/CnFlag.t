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


# for testing CnFlag class

use File::Copy qw(copy move);
use Text::Diff;

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::BrassFilter::CnFlag;

use Test::More 'no_plan';

use FindBin qw($Bin);

# existing entry
my $infile = $Bin.'/../testData/' . 'CnFlag_test.in';
my $infile_oeh = $Bin.'/../testData/' . 'CnFlag_test_oeh.in';
my $infile_ascat = $Bin.'/../testData/' . 'CnFlag_test_ascat.csv';
my $infile_ngs = $Bin.'/../testData/' . 'CnFlag_test_ngs.csv';
my $infile_bb = $Bin.'/../testData/' . 'CnFlag_test_bb.csv';
my $outfile = $Bin.'/../testData/' . 'CnFlag_test.out';
my $outfile_oeh = $Bin.'/../testData/' . 'CnFlag_test_oeh.out';
my $testfile = 'CnFlag_test.in';
my $testfile_oeh = 'CnFlag_test_oeh.in';

my $field = 25;
my $within = 500;

# get the infile ready
copy $infile, $testfile;
copy $infile_oeh, $testfile_oeh;

# make a new RG object
my $CnFlag = new Sanger::CGP::BrassFilter::CnFlag(-infile      => $testfile,
				     -one_end_hit => 0,
				     -ascat       => $infile_ascat,
				     -ngs         => $infile_ngs,
				     -bb          => $infile_bb,
				     -field       => $field,
                                     -within      => $within );
# process file
$CnFlag->process();

# make a new RG object
my $CnFlag_oeh = new Sanger::CGP::BrassFilter::CnFlag(-infile      => $testfile_oeh,
				     -one_end_hit => 1,
				     -ascat       => $infile_ascat,
				     -ngs         => $infile_ngs,
				     -bb          => $infile_bb,
				     -field       => $field,
                                     -within      => $within );
# process file
$CnFlag_oeh->process();

# check the size of the outfile is the same or greater
my $infile_size = -s $testfile;
my $outfile_size = -s $outfile;
my $infile_size_oeh = -s $testfile_oeh;
my $outfile_size_oeh = -s $outfile_oeh;

my $diff = diff "$testfile", "$outfile";
my $diff_oeh = diff "$testfile_oeh", "$outfile_oeh";

ok(defined $CnFlag, 'object defined');

is (($CnFlag->infile()), $testfile , "get infile");
is (($CnFlag->one_end_hit()), 0 , "one_end_hit");
is (($CnFlag->ascat()), $infile_ascat , "ascat_file");
is (($CnFlag->ngs()), $infile_ngs , "ngs_file");
is (($CnFlag->bb()),$infile_bb  , "bb_file");
is (($CnFlag->within()), $within , "within");
is ($infile_size, $outfile_size , "Outfile size check");
ok(defined $diff, 'file has been modified');

ok(defined $CnFlag_oeh, 'oeh object defined');

is (($CnFlag_oeh->infile()), $testfile_oeh , "oeh get infile");
is (($CnFlag_oeh->one_end_hit()), 1 , "oeh one_end_hit");
is (($CnFlag_oeh->ascat()), $infile_ascat , "oeh ascat_file");
is (($CnFlag_oeh->ngs()), $infile_ngs , "oeh ngs_file");
is (($CnFlag_oeh->bb()),$infile_bb  , "oeh bb_file");
is (($CnFlag_oeh->within()), $within , "oeh within");
is ($infile_size_oeh, $outfile_size_oeh , "oeh Outfile size check");
ok(defined $diff_oeh, 'oeh file has been modified');

if ($infile_size  == $outfile_size) { unlink $testfile; }
if ($infile_size_oeh  == $outfile_size_oeh) { unlink $testfile_oeh; }
