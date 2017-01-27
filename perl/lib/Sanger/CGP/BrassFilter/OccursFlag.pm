package Sanger::CGP::BrassFilter::OccursFlag;

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


# Author las
#
=head1 NAME

OccursFlag

=head1 SYNOPSIS

use Sanger::CGP::BrassFilter::OccursFlag;

my $OccursFlag = new Sanger::CGP::BrassFilter::OccursFlag(-infile   => $testfile,
					     -Lfield   => $Lfield,
					     -Hfield   => $Hfield,
					     -bin_size => $bin_size,
					     -within   => $within);

$OccursFlag->process();

=head1 DESCRIPTION

Class that updates the Occurrences flags on a filtered BrassI bedpe file.

That is, the number of time this lower or higher breakpoint coordinate occurs within this sample dataset.

bin_size - how to split up the dataset to manage memory.

It puts the number of L/H occurrences into the specified numbered fields.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

use strict;
use File::Copy qw(move);

use Bio::Brass;
our $VERSION = Bio::Brass->VERSION;

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::BrassFilter::OccursFlag();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};


    bless $self,$class;

    # defaults
    $self->{debug} = 0;
    $self->{Lfield} = 23; # which field of the file to put the balanced translocations in
    $self->{Hfield} = 24; # which field of the file to put the balanced translocations in
    $self->{bin_size} = 100000; # how to split up the dataset to manage memory
    $self->{within} = 500; # how far away the coordinates can be for it to qualify

    if ($args{-infile})   { $self->infile($args{-infile}); }
    if ($args{-Lfield})   { $self->Lfield($args{-Lfield}); }
    if ($args{-Hfield})   { $self->Hfield($args{-Hfield}); }
    if ($args{-bin_size}) { $self->bin_size($args{-bin_size}); }
    if ($args{-within})   { $self->within($args{-within}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 infile

  Arg (1)    : infile name
  Example    : $infile = $object->infile($infile);
  Description: name of the filtered brassI marked groups bedpe infile
  Return     : infile

=cut

sub infile {
    my $self = shift;
    $self->{infile} = shift if @_;
    return $self->{infile};
}
#-----------------------------------------------------------------------#

=head2 Lfield

  Arg (1)    : Lfield
  Example    : $Lfield = $object->field(field);
  Description: field of tab delimited file that the Loccurrences flag should be put into. default=23
  Return     : field

=cut

sub Lfield {
    my $self = shift;
    $self->{Lfield} = shift if @_;
    return $self->{Lfield};
}
#-----------------------------------------------------------------------#

=head2 Hfield

  Arg (1)    : Hfield
  Example    : $Hfield = $object->Hfield(Hfield);
  Description: field of tab delimited file that the Hoccurrences flag should be put into. default=24
  Return     : Hfield

=cut

sub Hfield {
    my $self = shift;
    $self->{Hfield} = shift if @_;
    return $self->{Hfield};
}
#-----------------------------------------------------------------------#

=head2 bin_size

  Arg (1)    : bin_size
  Example    : $bin_size = $object->bin_size(bin_size);
  Description: size of data bins. default=100000
  Return     : bin_size

=cut

sub bin_size {
    my $self = shift;
    $self->{bin_size} = shift if @_;
    return $self->{bin_size};
}
#-----------------------------------------------------------------------#

=head2 within

  Arg (1)    : within
  Example    : $within = $object->within(within);
  Description: how far away two coordinates can be for it to still be considered a hit. default=500
  Return     : within

=cut

sub within {
    my $self = shift;
    $self->{within} = shift if @_;
    return $self->{within};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 process

  Arg (0)    :
  Example    : $object->process();
  Description: process the infile and put the output in the field number specified
  Return     :

=cut

sub process {
    my ($self) = @_;

    # check the file is there and valid
    my $ok = $self->_check_file();
    unless ($ok) { print "Occurs: Check failed\n"; return; }

    $ok = $self->_read_data();
    unless ($ok) { print "Occurs: Read data failed\n"; return; }

    $ok = $self->_get_counts();
    unless ($ok) { print "Occurs: Get counts failed\n"; return; }

    $self->_print_file();
}
#-----------------------------------------------------------------------#
sub _check_file {
    my ($self) = @_;

    my $file = $self->{infile};

    unless ($file && (-e $file)) {
	print "file $file not found\n";
	return(0);
    }

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);

	my ($ok) = $self->_check_line($line);
	last;
    }
    close $fh;

    return(1);
}
#-----------------------------------------------------------------------#

sub _read_data {
    my $self = shift;

    my $file = $self->{infile};

    # get chr and L5/L3, H5/H3 into a hash with them binned (+/- 500 - edge ones go in 2 adjacent bins)
    # read all the rearrangements into memory
    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);
	chomp $line;

	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = $self->_check_line($line);
	return(0) unless ($chr1);

	$self->{rgs}->{$name}->{chr1} = $chr1;
	$self->{rgs}->{$name}->{start1} = $start1;
	$self->{rgs}->{$name}->{end1} = $end1;
	$self->_bin_occ_coords($chr1, $start1, $end1);
	$self->{rgs}->{$name}->{chr2} = $chr2;
	$self->{rgs}->{$name}->{start2} = $start2;
	$self->{rgs}->{$name}->{end2} = $end2;
	$self->_bin_occ_coords($chr2, $start2, $end2);
    }
    return(1);
}
#-----------------------------------------------------------------------#
sub _check_line {
    my ($self, $line) = @_;

    my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2);

    if ($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([\+-]?1?)\t([\+-]?1?)/) {
	($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
	return($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2);
    }
    else {
	print "entry: $line\nNot in bedpe format (chr1<TAB>start1<TAB>end1<TAB>chr2<TAB>start2<TAB>end2<TAB>name<TAB>score<TAB>strand1<TAB>strand2)\n";
	return(0);
    }
}
#-----------------------------------------------------------------------#
sub _bin_occ_coords {

    my ($self, $chr, $coord5, $coord3) = @_;

    my $bin_size = $self->{bin_size};
    my $within = $self->{within};

    my $div5 = int($coord5 / $bin_size);
    my $mod5 = $coord5 % $bin_size;
    my $bin5 = $div5 * $bin_size;

    my $div3 = int($coord3 / $bin_size);
    my $mod3 = $coord3 % $bin_size;
    my $bin3 = $div3 * $bin_size;

    my $key = $coord5 . '-' . $coord3;
    $self->{binned_coords}->{$chr}->{$bin5}->{$key}->{coord5} = $coord5;
    $self->{binned_coords}->{$chr}->{$bin5}->{$key}->{coord3} = $coord3;
    $self->{binned_coords}->{$chr}->{$bin3}->{$key}->{coord5} = $coord5;
    $self->{binned_coords}->{$chr}->{$bin3}->{$key}->{coord3} = $coord3;

    # if it is close to the lower coordinate end, add it to the preceeding bin too
    if ($mod5 <= ($within + 1)) {
	$self->{binned_coords}->{$chr}->{$bin3 - $bin_size}->{$key}->{coord5} = $coord5;
	$self->{binned_coords}->{$chr}->{$bin3 - $bin_size}->{$key}->{coord3} = $coord3;
    }
    # if it is close to the higher coordinate end, add it to the following bin too
    elsif (($bin_size - $mod3) <= ($within + 1)) {
	$self->{binned_coords}->{$chr}->{$bin5 + $bin_size}->{$key}->{coord5} = $coord5;
	$self->{binned_coords}->{$chr}->{$bin5 + $bin_size}->{$key}->{coord3} = $coord3;
    }
}
#-----------------------------------------------------------------------#

sub _get_counts {
    my $self = shift;

    # check through all rgs and find others in the set with similar coords
    foreach my $name(sort {$a <=> $b} keys %{$self->{rgs}}) {
	my $chr1 = $self->{rgs}->{$name}->{chr1};
	my $start1 = $self->{rgs}->{$name}->{start1};
	my $end1 = $self->{rgs}->{$name}->{end1};
	my $chr2 = $self->{rgs}->{$name}->{chr2};
	my $start2 = $self->{rgs}->{$name}->{start2};
	my $end2 = $self->{rgs}->{$name}->{end2};

	# get what bin L coord range in and
	# go through the contents of that bin and count how many things are within $within bases
	my $Lcount = $self->_get_occ_bin($chr1, $start1, $end1);
	$self->{rgs}->{$name}->{Lcount} = $Lcount;

	# get what bin H coord range in and
	# go through the contents of that bin and count how many things are within $within bases
	my $Hcount = $self->_get_occ_bin($chr2, $start2, $end2);
	$self->{rgs}->{$name}->{Hcount} = $Hcount;
    }
    return(1);
}
#-----------------------------------------------------------------------#
sub _get_occ_bin {
    my ($self, $chr, $coord5, $coord3) = @_;

    my $bin_size = $self->{bin_size};
    my $within = $self->{within};

    # get the coords bin
    my $bin = (int($coord5 / $bin_size)) * $bin_size;

    # go through that bin and count similar entries
    my $count = 0;
    foreach my $key(sort {$a cmp $b} keys %{$self->{binned_coords}->{$chr}->{$bin}}) {
	my $test5 = $self->{binned_coords}->{$chr}->{$bin}->{$key}->{coord5};
	my $test3 = $self->{binned_coords}->{$chr}->{$bin}->{$key}->{coord3};

	if ( ($coord5 < ($test3 + $within))  &&
             ($coord3 > ($test5 - $within)) ) {
	    $count++;
	}
    }

    return($count);
}
#-----------------------------------------------------------------------#

sub _print_file {
    my ($self) = @_;

    my $infile = $self->{infile};
    my $temp_file = $self->{infile} . '.temp';

    my $Lfield = $self->{Lfield} - 1; # arrays are zero base referenced
    my $Hfield = $self->{Hfield} - 1; # arrays are zero base referenced

    open my $fh, "<$infile" or die $!;
    open my $fh_temp, ">$temp_file" or die $!;

    while (my $line = <$fh>) {
	if ($line =~ /^\s*#/) { print $fh_temp $line; next; }
	next unless ($line);
	chomp $line;
	my @line = split "\t", $line;

	my $name = $line[6];

	if ($self->{rgs}->{$name}->{Lcount}) { $line[$Lfield] = $self->{rgs}->{$name}->{Lcount}; }
	if ($self->{rgs}->{$name}->{Hcount}) { $line[$Hfield] = $self->{rgs}->{$name}->{Hcount}; }
	my $done_line = join "\t", @line;

	print $fh_temp "$done_line\n";
    }
    close $fh;
    close $fh_temp;

    # check the size of the outfile is the same or greater
    my $infile_size = -s $infile;
    my $outfile_size = -s $temp_file;

    # move the new file to the old file name if the file is the expected size
    if ($outfile_size >= $infile_size) {
	move $temp_file, $infile;
    }
    else { print "WARN: translocation flagging failed!\n"; }
}

#-----------------------------------------------------------------------#

1;
