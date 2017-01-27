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


# for testing BlatFlag class

use File::Copy qw(copy move);
use Text::Diff;

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::BrassFilter::BlatFlag;
use File::Which qw(which);
use Test::More 'no_plan';
use FindBin qw($Bin);

# existing entry
my $infile = $Bin.'/../testData/' . 'BlatFlag_test.in';
my $outfile = $Bin.'/../testData/' . 'BlatFlag_test.out';
my $testfile = 'BlatFlag_test.in';
my $ref = '/nfs/cancer_ref01/human/37/genome.fa';
my $field = 26;
my $blat_script = which('blat');
my $minIdentity = 95;

# get the infile ready
copy $infile, $testfile;

# make a new RG object
my $BlatFlag = new Sanger::CGP::BrassFilter::BlatFlag(-infile => $testfile,
					 -field  => $field,
					 -blat   => $blat_script,
					 -minIdentity => $minIdentity,
					 -ref    => $ref );
# process file
$BlatFlag->process();

# check the size of the outfile is the same or greater
my $infile_size = -s $testfile;
my $outfile_size = -s $outfile;

my $diff = diff "$testfile", "$outfile";

ok(defined $BlatFlag, 'object defined');

is (($BlatFlag->infile()), $testfile , "get infile");
is (($BlatFlag->field()), $field , "field");
is (($BlatFlag->ref()), $ref , "ref");
ok(defined $diff, 'file has been modified');
if(-e '/nfs/cancer_ref01/human/37/genome.fa') {
  is ($infile_size, $outfile_size , "Outfile size check");
  if ($infile_size  == $outfile_size) { unlink $testfile; }
}
else {
  unlink $testfile;
}
