package Sanger::CGP::BrassFilter::CnFlag;

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


# Author las
#
=head1 NAME

CnFlag

=head1 SYNOPSIS

use Sanger::CGP::BrassFilter::CnFlag;

my $CnFlag = new Sanger::CGP::BrassFilter::CnFlag(-infile      => $testfile,
				     -one_end_hit => 0,
				     -ascat       => $infile_ascat,
				     -ngs         => $infile_ngs,
				     -bb          => $infile_bb,
				     -field       => $field,
                                     -within      => $within );

$CnFlag->process();


=head1 DESCRIPTION

Class that updates the near_copynumber_changepoint flags on a bedpe file.

3 different copy number segment files can be supplied.

It puts the flag string in the specified number field.

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
  Example    : $object = new Sanger::CGP::BrassFilter::CnFlag();
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
    $self->{field} = 25; # which field of the file to put the changepoints in
    $self->{within} = 100000; # how far away the coordinates can be for it to qualify

    if ($args{-infile})      { $self->infile($args{-infile}); }
    if ($args{-field})       { $self->field($args{-field}); }
    if ($args{-within})      { $self->within($args{-within}); }
    if ($args{-ascat})       { $self->ascat($args{-ascat}); }
    if ($args{-ngs})         { $self->ngs($args{-ngs}); }
    if ($args{-bb})          { $self->bb($args{-bb}); }
    if (defined($args{-one_end_hit})) { $self->one_end_hit($args{-one_end_hit}); }

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

=head2 field

  Arg (1)    : field
  Example    : $field = $object->field(field);
  Description: field of tab delimited file that the flag should be put into. default=25
  Return     : field

=cut

sub field {
    my $self = shift;
    $self->{field} = shift if @_;
    return $self->{field};
}
#-----------------------------------------------------------------------#

=head2 within

  Arg (1)    : within
  Example    : $within = $object->within(within);
  Description: how far away two coordinates can be for it to still be considered a hit. default=100000
  Return     : within

=cut

sub within {
    my $self = shift;
    $self->{within} = shift if @_;
    return $self->{within};
}
#-----------------------------------------------------------------------#

=head2 ascat

  Arg (1)    : $ascat (optional)
  Example    : $ascat = $pair->ascat($ascat);
  Description: ascat copynumber segments file
  Return     : $ascat

=cut

sub ascat {
    my $self = shift;
    $self->{ascat} = shift if @_;
    return $self->{ascat};
}
#-----------------------------------------------------------------------#

=head2 ngs

  Arg (1)    : $ngs (optional)
  Example    : $ngs = $pair->ngs($ngs);
  Description: ngs copynumber segments file
  Return     : $ngs

=cut

sub ngs {
    my $self = shift;
    $self->{ngs} = shift if @_;
    return $self->{ngs};
}
#-----------------------------------------------------------------------#

=head2 bb

  Arg (1)    : $bb (optional)
  Example    : $bb = $pair->bb($bb);
  Description: bb (Battenberg) copynumber segments file
  Return     : $bb

=cut

sub bb {
    my $self = shift;
    $self->{bb} = shift if @_;
    return $self->{bb};
}
#-----------------------------------------------------------------------#

=head2 one_end_hit

  Arg (1)    : 0/1
  Example    : $one_end_hit = $pair->one_end_hit($one_end_hit);
  Description: If set to 1, call a hit if only one end of the rearrangement is near a changepoint. default=0.
               Default behaviour (one_end_hit = 0) is to only call a hit if both ends are near a changepoint (of any sort, Affy, NGS or Battenberg)
  Return     : 0/1

=cut

sub one_end_hit {
    my $self = shift;
    $self->{one_end_hit} = shift if @_;
    return $self->{one_end_hit};
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
    my $ok = $self->_check_file($self->{infile});
    unless ($ok) { print "Cn: Check failed\n"; return; }

    # check optional copynumber files
    if ($self->{ascat}) {
	$ok = $self->_check_file($self->{ascat}, 'CN');
	unless ($ok) { print "Cn: Check ascat segments file failed\n"; return; }
    }
    if ($self->{ngs}) {
	$ok = $self->_check_file($self->{ngs}, 'CN');
	unless ($ok) { print "Cn: Check ngs segments file failed\n"; return; }
    }
    if ($self->{bb}) {
	$ok = $self->_check_file($self->{bb}, 'CN');
	unless ($ok) { print "Cn: Check battenberg segments file failed\n"; return; }
    }

    $ok = $self->_read_data();
    unless ($ok) { print "Cn: Read data failed\n"; return; }

    $ok = $self->_get_hits();
    unless ($ok) { print "Cn: Get hits failed\n"; return; }

    $self->_print_file();
}
#-----------------------------------------------------------------------#
sub _check_file {
    my ($self, $file, $cn) = @_;

    unless ($file && (-e "$file")) {
	print "file $file not found\n";
	return(0);
    }

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line =~ /\S/);

	if ($cn) {  $self->_check_cn_line($line); }
	else     {  $self->_check_line($line); }
	last;
    }
    close $fh;

    return(1);
}
#-----------------------------------------------------------------------#
sub _check_line {
    my ($self, $line) = @_;

    chomp $line;

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
sub _check_cn_line {
    my ($self, $line) = @_;

    chomp $line;

    my ($name,$chr,$start,$end,$Ntotal,$Nminor,$Ttotal,$Tminor);

    if ($line =~ /^.+?,\s*(\S+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*$/) {
	($chr,$start,$end,$Ntotal,$Nminor,$Ttotal,$Tminor) = ($1,$2,$3,$4,$5,$6,$7);
	return($chr,$start,$end,$Ntotal,$Nminor,$Ttotal,$Tminor);
    }
    elsif ($line =~ /^.+?,\s*(\S+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*$/) {
	($chr,$start,$end,$Ttotal,$Tminor) = ($1,$2,$3,$4,$5);
	return($chr,$start,$end,0,0,undef,undef);
    }
    else {
	print "entry: $line\nNot in copynumber segment format (name,chr,start,end,normal_total_copynumber(optional),normal_minor_copynumber(optional),tumour_total_copynumber,tumour_minor_copynumber)\n";
	return(0);
    }
}
#-----------------------------------------------------------------------#

sub _read_data {
    my $self = shift;

    if ($self->{ascat}) { $self->_get_breakpoints($self->{ascat},'ascat'); }
    if ($self->{ngs})  { $self->_get_breakpoints($self->{ngs},'ngs'); }
    if ($self->{bb})   { $self->_get_breakpoints($self->{bb},'bb'); }

    return(1);
}
#--------------------------------------------------------------------------------------#
sub _get_breakpoints {
    my ($self, $file, $type) = @_;

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line =~ /\S/);

	my ($chr,$start,$end,$Ntotal,$Nminor,$Ttotal,$Tminor) = $self->_check_cn_line($line);

	# if normal data exists (ie )
	# see if Normal total copynumber = Tumour total copynumber and Normal minor copynumber = Tumour minor copynumber
	# skip if they are the same (not interested in breakpoint that occur in both tumour and normal)
	next if (defined($Ttotal) && defined($Tminor) && ($Ntotal == $Ttotal) && ($Nminor == $Tminor));

	$self->{breakpoints}->{$chr}->{$start}->{$type} = 1;
	$self->{breakpoints}->{$chr}->{$end}->{$type} = 1;
    }
}
#--------------------------------------------------------------------------------------#
sub _get_hits {
    my $self = shift;

    $self->{hits} = {};

    my $file = $self->{infile};

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);

	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = $self->_check_line($line);
	return(0) unless ($chr1);
	my $L_does = $self->_match_to_change_point($chr1,$start1);
	my $H_does = $self->_match_to_change_point($chr2,$end2);

	# make sure there's only 1 result for each copynumber type is set for each end of the rearrangement
	my $L_merged_result_string = '';
	my $H_merged_result_string = '';
	if ($L_does =~ /A/)    { $L_merged_result_string .= 'A'; }
	if ($H_does =~ /A/)    { $H_merged_result_string .= 'A'; }
	if ($L_does =~ /N/)    { $L_merged_result_string .= 'N'; }
	if ($H_does =~ /N/)    { $H_merged_result_string .= 'N'; }
	if ($L_does =~ /B/)    { $L_merged_result_string .= 'B'; }
	if ($H_does =~ /B/)    { $H_merged_result_string .= 'B'; }
	my $merged_result_string = $L_merged_result_string . '_' . $H_merged_result_string;

	if ( ($self->{one_end_hit}) &&
	     ($L_does || $H_does) ) {
	    $self->{hits}->{$name} = $merged_result_string;
	}
	elsif ($L_does && $H_does) {
	    $self->{hits}->{$name} = $merged_result_string;
	}
	else {
	    $self->{hits}->{$name} = '_';
	}
    }
    close($fh);
    $self->{breakpoints} = undef;

    return(1);
}
#-----------------------------------------------------------------------------#
# internal method to check whether a position is near a changepoint (stored in a hash)
sub _match_to_change_point {
    my ($self, $chr, $pos) = @_;

    my $within = $self->{within};

    my $match = '';

    foreach my $change_point(keys %{$self->{breakpoints}->{$chr}}) {
	next unless (($pos > ($change_point - $within)) &&
                     ($pos < ($change_point + $within)));
	if ($self->{breakpoints}->{$chr}->{$change_point}->{ascat}){ $match .= 'A'; }
	if ($self->{breakpoints}->{$chr}->{$change_point}->{ngs})  { $match .= 'N'; }
	if ($self->{breakpoints}->{$chr}->{$change_point}->{bb})   { $match .= 'B'; }
    }
    return($match);
}

#-----------------------------------------------------------------------#

sub _print_file {
    my ($self) = @_;

    my $infile = $self->{infile};
    my $temp_file = $self->{infile} . '.temp';

    my $field = $self->{field} - 1; # arrays are zero base referenced

    open my $fh, "<$infile" or die $!;
    open my $fh_temp, ">$temp_file" or die $!;

    while (my $line = <$fh>) {
	if ($line =~ /^\s*#/) { print $fh_temp $line; next; }
	next unless ($line);
	chomp $line;
	my @line = split "\t", $line;

	my $name = $line[6];

	$line[$field] = $self->{hits}->{$name};

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
    else { print "WARN: Copynumber changepoint flagging failed!\n"; }
}

#-----------------------------------------------------------------------#

1;
