package Sanger::CGP::BrassFilter::BlatFlag;

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


#
# Author las
#
=head1 NAME

BlatFlag

=head1 SYNOPSIS

use Sanger::CGP::BrassFilter::BlatFlag;

my $BlatFlag = new Sanger::CGP::BrassFilter::BlatFlag(-infile => $testfile,
					 -field  => $field,
					 -blat   => $blat_script,
					 -minIdentity => $minIdentity,
					 -ref    => $ref );
# process file
$BlatFlag->process();

=head1 DESCRIPTION

Class that updates the Lrange to Hrange blat score flag on a bedpe file. It puts the score in the specified number field.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

use strict;
use File::Copy qw(move);
use Bio::DB::HTS;
use File::Temp qw(tempdir);
use File::Spec;

use Bio::Brass;
our $VERSION = Bio::Brass->VERSION;

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::BrassFilter::BlatFlag();
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
    $self->{field} = 26; # which field of the file to put the blat score in
    $self->{blat} = 'blat';
    $self->{minIdentity} = 95;

    if ($args{-infile})      { $self->infile($args{-infile}); }
    if ($args{-field})       { $self->field($args{-field}); }
    if ($args{-ref})         { $self->ref($args{-ref}); }
    if ($args{-blat})        { $self->blat($args{-blat}); }
    if ($args{-minIdentity}) { $self->minIdentity($args{-minIdentity}); }

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
  Description: field of tab delimited file that the flag should be put into. default=26
  Return     : field

=cut

sub field {
    my $self = shift;
    $self->{field} = shift if @_;
    return $self->{field};
}
#-----------------------------------------------------------------------#

=head2 ref

  Arg (1)    : ref
  Example    : $ref = $object->ref(ref);
  Description: reference genome in fasta format, required to retrieve L and H range sequence.
               fai index must also be present.
  Return     : ref

=cut

sub ref {
    my $self = shift;
    $self->{ref} = shift if @_;
    return $self->{ref};
}
#-----------------------------------------------------------------------#

=head2 blat

  Arg (1)    : blat
  Example    : $blat = $object->blat(blat);
  Description: blat script (include path if necessary). default=blat
  Return     : blat

=cut

sub blat {
    my $self = shift;
    $self->{blat} = shift if @_;
    return $self->{blat};
}
#-----------------------------------------------------------------------#

=head2 minIdentity

  Arg (1)    : minIdentity
  Example    : $minIdentity = $object->minIdentity(minIdentity);
  Description: minIdentity value to supply to blat. default=95
  Return     : minIdentity

=cut

sub minIdentity {
    my $self = shift;
    $self->{minIdentity} = shift if @_;
    return $self->{minIdentity};
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
    unless ($ok) { print "Blat: Check failed\n"; return; }

    # check the ref file is there and valid
    $ok = $self->_check_file($self->{ref}, 'ref');
    unless ($ok) { print "Blat: Check ref failed\n"; return; }

    $ok = $self->_read_data();
    unless ($ok) { print "Blat: Read data failed\n"; return; }

    $ok = $self->_get_hits();
    unless ($ok) { print "Blat: Get hits failed\n"; return; }

    $self->_print_file();
}
#-----------------------------------------------------------------------#
sub _check_file {
    my ($self, $file, $ref) = @_;

    unless ($file && (-e $file)) {
	print "file $file not found\n";
	return(0);
    }

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line =~ /\S/);

	if ($ref) { $self->_check_ref_line($line); }
	else      { $self->_check_line($line); }
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
sub _check_ref_line {
    my ($self, $line) = @_;

    chomp $line;

    return($line);
}
#-----------------------------------------------------------------------#

sub _read_data {
    my $self = shift;

    my $file = $self->{infile};

    my $data = {};

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);

	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = $self->_check_line($line);
	return(0) unless ($chr1);

	$self->{data}->{$name}->{Lrange} = "$chr1:" . ($start1 + 1) . "-$end1";
	$self->{data}->{$name}->{Hrange} = "$chr2:" . ($start2 + 1) . "-$end2";
    }
    close($fh);

    return(1);
}
#--------------------------------------------------------------------------------------#

sub _get_hits {
  my $self = shift;

  my $blat = $self->{blat};
  my $minIdentity = $self->{minIdentity};

  my $tempdir = tempdir( 'BlatFlagXXXXXX', DIR => File::Spec->tmpdir(), CLEANUP => 1 );
  my $Lfile = "$tempdir/Lrange.fa";
  my $Hfile = "$tempdir/Hrange.fa";

  # load an indexed fasta file
  my $fai = Bio::DB::HTS::Fai->load( $self->{ref} );

  foreach my $name(keys %{$self->{data}}) {
    my $Lrange = $self->{data}->{$name}->{Lrange};
    my $Hrange = $self->{data}->{$name}->{Hrange};

    open my $fhl, ">$Lfile" or die $!;
    print $fhl sprintf ">%s\n%s\n", $Lrange, $fai->fetch($Lrange);
    close $fhl;

    open my $fhh, ">$Hfile" or die $!;
    print $fhh sprintf ">%s\n%s\n", $Hrange, $fai->fetch($Hrange);
    close $fhh;

    my $command = "$blat $Lfile $Hfile -noHead -minIdentity=$minIdentity /dev/stdout";
    my $score = 0;
    my ($pid, $process);
    $pid = open $process, q{-|}, $command or die 'Could not fork: '.$!;
    while(my $line = <$process>) {
      next unless($line =~ m/$Lrange/);
      my @hit = split /\t/, $line; # take the top blat hit
      $score = $hit[0] - $hit[1];
      last;
    }
    close $process;

    $self->{data}->{$name}->{score} = $score;
    if ($self->{debug}) { print "$name | SCORE:$score\n"; }
  }
  # save some file ops, only delete the L/H files at end of loop
	unlink $Lfile or die $! if(-e $Lfile);
	unlink $Hfile or die $! if(-e $Lfile);

  return(1);
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

	$line[$field] = $self->{data}->{$name}->{score};

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
    else { print "WARN: Blat flagging failed!\n"; }
}

#-----------------------------------------------------------------------#

1;
