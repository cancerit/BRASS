package Sanger::CGP::BrassFilter::TransFlag;

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


# Author las
#
=head1 NAME

TransFlag

=head1 SYNOPSIS

use Sanger::CGP::BrassFilter::TransFlag;

my $TransFlag = new Sanger::CGP::BrassFilter::TransFlag(-infile => $testfile,
                                           -bal_trans_field => $bal_trans_field,
                                           -inv_field => $inv_field,
                                           -bal_distance => $bal_distance,
                                           -inv_distance => $inv_distance);

$TransFlag->process();

=head1 DESCRIPTION

Class that updates the translocation flags on a filtered BrassI bedpe file.

It puts the ids (comma separated string) of paired inversion/balanced_translocation breakpoints into the specified numbered fields.

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
  Example    : $object = new Sanger::CGP::BrassFilter::TransFlag();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};


    bless $self,$class;

    $self->{debug} = 0;
    $self->{bal_trans_field} = 21; # which field of the file to put the balanced translocations in
    $self->{inv_field} = 22; # which field of the file to put the inversions in
    $self->{bal_distance} = 100000; # how far away the coordinates can be for it to qualify as a balanced translocation
    $self->{inv_distance} = 1000; # how far away the coordinates can be for it to qualify as an inversion

    if ($args{-infile})          { $self->infile($args{-infile}); }
    if ($args{-bal_trans_field}) { $self->bal_trans_field($args{-bal_trans_field}); }
    if ($args{-inv_field})       { $self->inv_field($args{-inv_field}); }
    if ($args{-bal_distance})    { $self->bal_distance($args{-bal_distance}); }
    if ($args{-inv_distance})    { $self->inv_distance($args{-inv_distance}); }

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

=head2 bal_trans_field

  Arg (1)    : bal_trans_field
  Example    : $bal_trans_field = $object->field(field);
  Description: field of tab delimited file that the balanced translocation flag should be put into. default=21
  Return     : field

=cut

sub bal_trans_field {
    my $self = shift;
    $self->{bal_trans_field} = shift if @_;
    return $self->{bal_trans_field};
}
#-----------------------------------------------------------------------#

=head2 inv_field

  Arg (1)    : inv_field
  Example    : $inv_field = $object->inv_field(inv_field);
  Description: field of tab delimited file that the inversion flag should be put into. default=22
  Return     : inv_field

=cut

sub inv_field {
    my $self = shift;
    $self->{inv_field} = shift if @_;
    return $self->{inv_field};
}
#-----------------------------------------------------------------------#

=head2 bal_distance

  Arg (1)    : bal_distance
  Example    : $bal_distance = $object->bal_distance(bal_distance);
  Description: how far away two coordinates can be for it to still be considered a balanced translocation. default=100000
  Return     : bal_distance

=cut

sub bal_distance {
    my $self = shift;
    $self->{bal_distance} = shift if @_;
    return $self->{bal_distance};
}
#-----------------------------------------------------------------------#

=head2 inv_distance

  Arg (1)    : inv_distance
  Example    : $inv_distance = $object->inv_distance(inv_distance);
  Description: how far away two coordinates can be for it to still be considered an inversion. default=1000
  Return     : inv_distance

=cut

sub inv_distance {
    my $self = shift;
    $self->{inv_distance} = shift if @_;
    return $self->{inv_distance};
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
    return unless $ok;

    $ok = $self->_read_data();
    return unless $ok;

    $ok = $self->_get_bals();
    return unless $ok;

    $ok = $self->_get_invs();
    return unless $ok;

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
    my ($self) = @_;

    my $file = $self->{infile};

   # read all the rearrangements into memory
    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);
	chomp $line;

	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = $self->_check_line($line);
	return(0) unless ($chr1);

	if ($strand1 eq '-1')     { $strand1 = '-'; }
	elsif ($strand1 =~ /\+?1/) { $strand1 = '+'; }
	if ($strand2 eq '-1')     { $strand2 = '-'; }
	elsif ($strand2 =~ /\+?1/) { $strand2 = '+'; }

	# load into inversion data hash
	if ($chr1 eq $chr2) {
	    if (($strand1 eq '+') && ($strand2 eq '+')) {
		$self->{inv_data}->{$chr1}->{'F'}->{$name}->{start1} = $start1;
		$self->{inv_data}->{$chr1}->{'F'}->{$name}->{end1} = $end1;
		$self->{inv_data}->{$chr1}->{'F'}->{$name}->{start2} = $start2;
		$self->{inv_data}->{$chr1}->{'F'}->{$name}->{end2} = $end2;
	    }
	    elsif (($strand1 eq '-') && ($strand2 eq '-')) {
		$self->{inv_data}->{$chr1}->{'R'}->{$name}->{start1} = $start1;
		$self->{inv_data}->{$chr1}->{'R'}->{$name}->{end1} = $end1;
		$self->{inv_data}->{$chr1}->{'R'}->{$name}->{start2} = $start2;
		$self->{inv_data}->{$chr1}->{'R'}->{$name}->{end2} = $end2;
	    }
	}
	# load into balanced translocation data hash
	else {
	    my $chr_pair = $chr1 . '_' . $chr2;
	    $self->{bal_data}->{$chr_pair}->{$name}->{chr1} = $chr1;
	    $self->{bal_data}->{$chr_pair}->{$name}->{start1} = $start1;
	    $self->{bal_data}->{$chr_pair}->{$name}->{end1} = $end1;
	    $self->{bal_data}->{$chr_pair}->{$name}->{strand1} = $strand1;
	    $self->{bal_data}->{$chr_pair}->{$name}->{chr2} = $chr2;
	    $self->{bal_data}->{$chr_pair}->{$name}->{start2} = $start2;
	    $self->{bal_data}->{$chr_pair}->{$name}->{end2} = $end2;
	    $self->{bal_data}->{$chr_pair}->{$name}->{strand2} = $strand2;
	}

    }
    close $fh;

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

sub _get_bals {
    my $self = shift;

    my $data = $self->{bal_data}; # data from file in memory
    my $distance = $self->{bal_distance}; # how far away the coordinates can be


    foreach my $chr_pair(keys %{$self->{bal_data}}) {
	foreach my $name1(sort keys %{$self->{bal_data}->{$chr_pair}}) {
	    foreach my $name2(sort keys %{$self->{bal_data}->{$chr_pair}}) {
		next if ($name1 eq $name2);

		my $chr1_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{chr1};
		my $start1_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{start1};
		my $end1_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{end1};
		my $strand1_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{strand1};
		my $chr2_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{chr2};
		my $start2_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{start2};
		my $end2_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{end2};
		my $strand2_1 = $self->{bal_data}->{$chr_pair}->{$name1}->{strand2};

		my $chr1_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{chr1};
		my $start1_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{start1};
		my $end1_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{end1};
		my $strand1_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{strand1};
		my $chr2_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{chr2};
		my $start2_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{start2};
		my $end2_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{end2};
		my $strand2_2 = $self->{bal_data}->{$chr_pair}->{$name2}->{strand2};

		if (($chr1_1 eq $chr1_2) &&
                    ($chr2_1 eq $chr2_2) &&
                    ($strand1_1 ne $strand1_2) &&
                    ($strand2_1 ne $strand2_2) &&
                    ($start1_2 > ($start1_1 - $distance)) &&
                    ($end1_2   < ($end1_1   + $distance)) &&
                    ($start2_2 > ($start2_1 - $distance)) &&
                    ($end2_2   < ($end2_1   + $distance)) ) {

		    if ($self->{bal_data}->{$chr_pair}->{$name1}->{res}) { $self->{bal_data}->{$chr_pair}->{$name1}->{res} .= ",$name2"; }
		    else                                                 { $self->{bal_data}->{$chr_pair}->{$name1}->{res} = $name2; }
		}
                # if swapped over...
		elsif (($chr1_1 eq $chr2_2) &&
                       ($chr2_1 eq $chr1_2) &&
                       ($strand1_1 ne $strand2_2) &&
                       ($strand2_1 ne $strand1_2) &&
                       ($start1_2 > ($start2_1 - $distance)) &&
                       ($end1_2   < ($end2_1   + $distance)) &&
                       ($start2_2 > ($start1_1 - $distance)) &&
                       ($end2_2   < ($end1_1   + $distance)) ) {
		    if ($self->{bal_data}->{$chr_pair}->{$name1}->{res}) { $self->{bal_data}->{$chr_pair}->{$name1}->{res} .= ",$name2"; }
		    else                                                 { $self->{bal_data}->{$chr_pair}->{$name1}->{res} = $name2; }
		}

	    }
	}
    }

    return(1);

    # look for pairs of rearrangements with the opposite strands but same coordinates
    # (take into account that L and H might be swapped over)
    # takes about 30 sec
    my $query = "select distinct rg1.id_rearrangement_group
                                ,rg2.id_rearrangement_group
                 from rg_rearrangement_group rg1
                     ,rg_rearrangement_group rg2
                 where rg1.id_analysis_proc = ?
                 and rg2.id_analysis_proc = ?
                 and rg1.paired_flag = 32
                 and nvl(rg1.repeat, -1) != 'newDR'
                 and rg2.paired_flag = 32
                 and  nvl(rg2.repeat, -1) != 'newDR'

                 and ( ( rg1.chrL = rg2.chrL
                         and rg1.chrH = rg2.chrH
                         and rg1.strandL != rg2.strandL
                         and rg1.strandH != rg2.strandH
                         and rg2.L5 > (rg1.L5 - $distance)
                         and rg2.L3 < (rg1.L3 + $distance)
                         and rg2.H5 > (rg1.H5 - $distance)
                         and rg2.H3 < (rg1.H3 + $distance)
                       )

                       or
                       (  rg1.chrL = rg2.chrH
                         and rg1.chrH = rg2.chrL
                         and rg1.strandL != rg2.strandH
                         and rg1.strandH != rg2.strandL
                         and rg2.L5 > (rg1.H5 - $distance)
                         and rg2.L3 < (rg1.H3 + $distance)
                         and rg2.H5 > (rg1.L5 - $distance)
                         and rg2.H3 < (rg1.L3 + $distance)
                       )
                     )
                  order by rg1.id_rearrangement_group, rg2.id_rearrangement_group";

}
#-----------------------------------------------------------------------#

sub _get_invs {
    my $self = shift;

    my $data = $self->{inv_data}; # dataset in memory
    # by default, look for read pair 2 (type 8) within 2 insert sizes of the inversion size
    my $distance = $self->{inv_distance}; # how far away the coordinates can be


    foreach my $chr(keys %{$self->{inv_data}}) {
	foreach my $nameF(sort keys %{$self->{inv_data}->{$chr}->{'F'}}) {
	    foreach my $nameR(sort keys %{$self->{inv_data}->{$chr}->{'R'}}) {
		my $start1F = $self->{inv_data}->{$chr}->{'F'}->{$nameF}->{start1};
		my $start1R = $self->{inv_data}->{$chr}->{'R'}->{$nameR}->{start1};
		my $start2F = $self->{inv_data}->{$chr}->{'F'}->{$nameF}->{start2};
		my $start2R = $self->{inv_data}->{$chr}->{'R'}->{$nameR}->{start2};

		if ( ( (($start1R + $distance) > $start1F) &&
                       (($start2R + $distance) > $start2F) &&
                       (($start1F + $distance) > $start1R) &&
                       (($start2F + $distance) > $start2R) ) ||
                     ( (($start1R + $distance) > $start2F) &&
                       (($start2R + $distance) > $start1F) &&
                       (($start1F + $distance) > $start2R) &&
                       (($start2F + $distance) > $start1R) ) ) {

#		print "INV: $chr $nameF $nameR $start1F $start1R $start2F $start2R\n";

		    if ($self->{inv_data}->{$chr}->{'F'}->{$nameF}->{res}) { $self->{inv_data}->{$chr}->{'F'}->{$nameF}->{res} .= ",$nameR"; }
		    else                                                   { $self->{inv_data}->{$chr}->{'F'}->{$nameF}->{res} = $nameR; }
		    if ($self->{inv_data}->{$chr}->{'R'}->{$nameR}->{res}) { $self->{inv_data}->{$chr}->{'R'}->{$nameR}->{res} .= ",$nameF"; }
		    else                                                   { $self->{inv_data}->{$chr}->{'R'}->{$nameR}->{res} = $nameF; }

		}
	    }
	}
    }
    return(1);
}
#-----------------------------------------------------------------------#

sub _print_file {
    my ($self) = @_;

    my $infile = $self->{infile};
    my $temp_file = $self->{infile} . '.temp';

    my $bal_field = $self->{bal_trans_field} - 1; # arrays are zero base referenced
    my $inv_field = $self->{inv_field} - 1; # arrays are zero base referenced

    open my $fh, "<$infile" or die $!;
    open my $fh_temp, ">$temp_file" or die $!;

    while (my $line = <$fh>) {
	if ($line =~ /^\s*#/) { print $fh_temp $line; next; }
	next unless ($line);
	chomp $line;
	my @line = split "\t", $line;

	my $name = $line[6];
	my $chr1 = $line[0];
	my $chr2 = $line[3];
	my $chr_pair1 = $chr1 . '_' . $chr2;
	my $chr_pair2 = $chr2 . '_' . $chr1;

	if ($self->{inv_data}->{$chr1}->{'F'}->{$name}->{res}) { $line[$inv_field] = $self->{inv_data}->{$chr1}->{'F'}->{$name}->{res}; }
	if ($self->{inv_data}->{$chr1}->{'R'}->{$name}->{res}) { $line[$inv_field] = $self->{inv_data}->{$chr1}->{'R'}->{$name}->{res}; }


	if ($self->{bal_data}->{$chr_pair1}->{$name}->{res} &&
            $self->{bal_data}->{$chr_pair2}->{$name}->{res}) {
	    $line[$bal_field] = $self->{bal_data}->{$chr_pair1}->{$name}->{res} . ',' . $self->{bal_data}->{$chr_pair2}->{$name}->{res};
	}
	elsif ($self->{bal_data}->{$chr_pair1}->{$name}->{res}) { $line[$bal_field] = $self->{bal_data}->{$chr_pair1}->{$name}->{res}; }
	elsif ($self->{bal_data}->{$chr_pair2}->{$name}->{res}) { $line[$bal_field] = $self->{bal_data}->{$chr_pair2}->{$name}->{res}; }

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
