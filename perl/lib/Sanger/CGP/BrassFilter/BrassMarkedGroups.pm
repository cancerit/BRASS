package Sanger::CGP::BrassFilter::BrassMarkedGroups;

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

use Bio::Brass;
our $VERSION = Bio::Brass->VERSION;

# Author las
#
=head1 NAME

BrassMarkedGroups

=head1 SYNOPSIS

use Sanger::CGP::BrassFilter::BrassMarkedGroups;

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

=head1 DESCRIPTION

Class that handles filtering and outputting of a BrassI marked groups file. Output is in bedpe format.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::BrassFilter::BrassMarkedGroups();
  Description: make a new BrassMarkedGroup object
  Return     : SumFile object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};


    bless $self,$class;

    $self->{debug} = $args{-debug} || 0;

    # default filter values (if not passed in)
    $self->{min_tumour_count_high} = 4;
    $self->{min_tumour_count_low} = 2;
    $self->{max_normal_count} = 0;
    $self->{max_np_count} = 0;
    $self->{max_np_sample_count} = 0;
    $self->{distance_threshold} = 100;
    $self->{seq_depth_threshold} = 25;
    $self->{seq_depth} = 30;
    $self->{discard_if_repeats} = 0;

    if ($args{-infile})  { $self->infile($args{-infile}); }
    if ($args{-outfile}) { $self->outfile($args{-outfile}); }

    # filters
    if (defined($args{-tumour}))                { $self->tumour($args{-tumour}); }
    if (defined($args{-min_tumour_count_high})) { $self->min_tumour_count_high($args{-min_tumour_count_high}); }
    if (defined($args{-min_tumour_count_low}))  { $self->min_tumour_count_low($args{-min_tumour_count_low}); }
    if (defined($args{-max_normal_count}))     { $self->max_normal_count($args{-max_normal_count}); }
    if (defined($args{-max_np_count}))         { $self->max_np_count($args{-max_np_count}); }
    if (defined($args{-max_np_sample_count}))  { $self->max_np_sample_count($args{-max_np_sample_count}); }
    if (defined($args{-distance_threshold}))   { $self->distance_threshold($args{-distance_threshold}); }
    if (defined($args{-seq_depth_threshold}))  { $self->seq_depth_threshold($args{-seq_depth_threshold}); }
    if (defined($args{-seq_depth}))            { $self->seq_depth($args{-seq_depth}); }
    if (defined($args{-discard_if_repeats}))   { $self->discard_if_repeats($args{-discard_if_repeats}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 outfile

  Arg (0)    : $outfile
  Example    : $outfile = $object->outfile($outfile);
  Description: required outfile  (bedpe filename extension will be appended if not supplied)
  Return     : $outfile

=cut

sub outfile {
    my $self = shift;
    $self->{outfile} = shift if @_;
    return $self->{outfile};
}
#-----------------------------------------------------------------------#

=head2 infile

  Arg (0)    : infile name
  Example    : $infile = $object->infile($infile);
  Description: name of the brassI marked groups infile
  Return     : infile

=cut

sub infile {
    my $self = shift;
    my $infile = shift;
    if ($infile) {
	$self->{infile} = $infile;
	my $count = `wc -l $self->{infile}`;
	if ($count =~ /^\s*(\d+)\s+/) { $self->{file_count} = $1; }
    }
    return $self->{infile};
}
#-----------------------------------------------------------------------#

=head2 tumour

  Arg (0)    : $tumour
  Example    : $tumour = $object->tumour($tumour);
  Description: tumour sample name
  Return     : $tumour

=cut

sub tumour {
    my $self = shift;
    $self->{tumour} = shift if @_;
    return $self->{tumour};
}
#-----------------------------------------------------------------------#

=head2 file_count

  Arg (0)    : file_count
  Example    : $file_count = $object->file_count();
  Description: getter for file_count - the number of rows in the file
  Return     : file_count

=cut

sub file_count {
    my $self = shift;
    return $self->{file_count};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 header

  Arg (0)    : header
  Example    : $header = $object->header();
  Description: file header from brass I input file
  Return     : header

=cut

sub header {
    my $self = shift;
    $self->{header} = shift if @_;
    return $self->{header};
}
#-----------------------------------------------------------------------#
#### filters ######
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 min_tumour_count_low

  Arg (1)    : $min_tumour_count_low
  Example    : $min_tumour_count_low = $object->min_tumour_count_low($min_tumour_count_low);
  Description: only load entries with at least this number of tumour sample reads.
               *This flag if applied if seq_depth_threshold is not exceeded*
  Return     : $min_tumour_count_low

=cut

sub min_tumour_count_low {
    my $self = shift;
    $self->{min_tumour_count_low} = shift if @_;
    return $self->{min_tumour_count_low};
}
#-----------------------------------------------------------------------#

=head2 min_tumour_count_high

  Arg (1)    : $min_tumour_count_high
  Example    : $min_tumour_count_high = $object->min_tumour_count_high($min_tumour_count_high);
  Description: only load entries with at least this number of tumour sample reads.
               *This flag is applied if seq_depth_threshold is exceeded*
  Return     : $min_tumour_count_high

=cut

sub min_tumour_count_high {
    my $self = shift;
    $self->{min_tumour_count_high} = shift if @_;
    return $self->{min_tumour_count_high};
}
#-----------------------------------------------------------------------#

=head2 max_normal_count

  Arg (1)    : $max_normal_count
  Example    : $max_normal_count = $object->max_normal_count($max_normal_count);
  Description: only load entries with at most this number of normal sample reads
  Return     : $max_normal_count

=cut

sub max_normal_count {
    my $self = shift;
    $self->{max_normal_count} = shift if @_;
    return $self->{max_normal_count};
}
#-----------------------------------------------------------------------#

=head2 max_np_count

  Arg (1)    : $max_np_count
  Example    : $max_np_count = $object->max_np_count($max_np_count);
  Description: only load entries with at most this number of normal panel sample reads
  Return     : $max_np_count

=cut

sub max_np_count {
    my $self = shift;
    $self->{max_np_count} = shift if @_;
    return $self->{max_np_count};
}
#-----------------------------------------------------------------------#

=head2 max_np_sample_count

  Arg (1)    : $max_np_sample_count
  Example    : $max_np_sample_count = $object->max_np_sample_count($max_np_sample_count);
  Description: only load entries where up to this number of np samples have reads
  Return     : $max_np_count

=cut

sub max_np_sample_count {
    my $self = shift;
    $self->{max_np_sample_count} = shift if @_;
    return $self->{max_np_sample_count};
}
#-----------------------------------------------------------------------#

=head2 seq_depth_threshold

  Arg (1)    : $seq_depth_threshold
  Example    : $seq_depth_threshold = $object->seq_depth_threshold($seq_depth_threshold);
  Description: sequence depth threshold - if the samples seq_depth is above this, min_tumour_count_high is applied in place of min_tumour_count_low
  Return     : $seq_depth_threshold

=cut

sub seq_depth_threshold {
    my $self = shift;
    $self->{seq_depth_threshold} = shift if @_;
    return $self->{seq_depth_threshold};
}
#-----------------------------------------------------------------------#

=head2 distance_threshold

  Arg (1)    : $distance_threshold
  Example    : $distance_threshold = $object->distance_threshold($distance_threshold);
  Description: distance between lower and upper breakpoints threshold - if breakpoints are on the same chromosome, only process them if over this value
  Return     : $distance_threshold

=cut

sub distance_threshold {
    my $self = shift;
    $self->{distance_threshold} = shift if @_;
    return $self->{distance_threshold};
}
#-----------------------------------------------------------------------#

=head2 seq_depth

  Arg (1)    : $seq_depth
  Example    : $seq_depth = $object->seq_depth($seq_depth);
  Description: sequence depth (default = 30x).
  Return     : $seq_depth

=cut

sub seq_depth {
    my $self = shift;
    $self->{seq_depth} = shift if @_;
    return $self->{seq_depth};
}
#-----------------------------------------------------------------------#

=head2 discard_if_repeats

  Arg (1)    : $discard_if_repeats
  Example    : $discard_if_repeats = $object->discard_if_repeats($discard_if_repeats);
  Description: dont load any entries with repeats
  Return     : $discard_if_repeats

=cut

sub discard_if_repeats {
    my $self = shift;
    $self->{discard_if_repeats} = shift if @_;
    return $self->{discard_if_repeats};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 samples

  Arg (0)    : samples
  Example    : $samples = $object->samples($samples);
  Description: returns a reference to an array of samples that are represented in this summary file
  Return     : samples

=cut

sub samples {
    my $self = shift;
    my $samples = shift if @_;

    #
    unless ($self->{samples}) { $self->get_sample_types(); }

    my @sample_list = ();
    foreach (@{$self->{samples}}) { push @sample_list, $_->[0]; }

    return \@sample_list;
}

#-----------------------------------------------------------------------#

=head2 sample_type

  Arg (0)    :
  Example    : $sample_type = $object->sample_type($sample);
  Description: returns a sample_type for a given sample (that is listed in samples)
  Return     : sample_type

=cut

sub sample_type {
    my ($self, $sample) = @_;
    my $sample_type = '';

    foreach my $entry(@{$self->{samples}}) {
	next unless ($entry->[0] eq $sample);
	$sample_type = $entry->[1];
	last;
    }

    return($sample_type);
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 get_sample_types

  Arg ()     :
  Example    : $object->get_sample_types();
  Description: gets the sample_type (T = tumour, N = normal, M = metastasis, P = normal panel) for a sample synonym from brassI marked groups file headers
               (Note the order of samples in the header is reversed in some versions of brass - the code handles this)
  Return     : $sample_types

=cut

sub get_sample_types {
    my $self = shift;

    my $file = $self->{infile};
    return unless ($file && -e $file);

    # get the #SAMPLE lines out of the header
    my $output = `grep #SAMPLE $file`;
    my @lines = split "\n", $output;

   # work out sample_types
    my $tumour_found = 0;
    my $normal_found = 0;
    foreach my $line(@lines) {
	my $sample = '';
	my $type = '';
	if ($line =~ /SAMPLE\s+\d+\s+(\S+)\s*$/) { $sample = $1; }

	if    ($sample =~ /^NP_/)     { $type = 'P'; }
	elsif ($self->{tumour} eq $sample) { $type = 'T'; $tumour_found = 1; }
	else { $type = 'N'; }

        push @{$self->{samples}}, [$sample, $type];
    }

    if ($self->{debug}) { foreach (@{$self->{samples}}) { print $_->[0] . q{ } . $_->[1] . "\n"; } }
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 process

  Arg ()     :
  Example    : $object->process();
  Description: runs the filtering process and outputs to output file
  Return     : $count of lines processed

=cut

sub process {
  my $self = shift;

  my $infile = $self->{infile};

  unless ($infile && -e $infile) {
    print "can not process. brassI marked groups infile (infile) not set or not found\n";
    die $!;
  }

  open my $fh, "<$infile" or die $!;
  my $chunk = 1;
  my $count = 0;

  $self->_parse_header($fh); # get the list of samples involved out of the header block

  # set the min_tumour count to use (depends on seq_depth_threshold)
  $self->{min_tumour_count} = $self->{min_tumour_count_low};

  if ($self->{seq_depth} > $self->{seq_depth_threshold}) {
    $self->{min_tumour_count} = $self->{min_tumour_count_high};
  }

  print "min_tumour_count $self->{min_tumour_count}\n" if($self->{debug});

  # Print header to file here rather than in print row.
  # That way we end up with an empty file with header so we know the script ran OK.

  # open an output file handler if it isn't already set
  my $bedpe_file = ($self->{outfile} || $self->{infile});
  unless ($bedpe_file =~ /\.bedpe$/) { $bedpe_file .= '.bedpe'; }
  my $bedpe_fh = $self->{bedpe_out_fh};

  unless ($bedpe_fh) {
    if ($self->{debug}) { print "opening $bedpe_file\n"; }
    open $bedpe_fh, "> $bedpe_file" or die $!;
    $self->{bedpe_out_fh} = $bedpe_fh;
    # print the header
    print $bedpe_fh $self->{header};
    print $bedpe_fh "# chr1\tstart1\tend1\tchr2\tstart2\tend2\tid/name\tbrass_score\tstrand1\tstrand2\trepeats\tnp_sample_count\ttumour_count\tnormal_count\tnp_count\tbkpt_distance\tsample\tsample_type\tnames\tcount\tbal_trans\tinv\toccL\toccH\tcopynumber_flag\trange_blat\n";
  }

  my $check_count = 0;
  while (my $line = <$fh>) {
    chomp $line;
    next unless ($line =~ /\S+\s+[+-]\s+\d+\s+\d+\s+\S+\s+[+-]\s+\d+\s+\d+\s+/); # chrL,strandL,L5,L3,chrH,strandH,H5,H3
    my @line = split "\t", $line;
    $count++ if($self->_add_row(\@line, ($count + 1)) == 1);
  }

  close $fh;
  close $self->{bedpe_out_fh} if ($self->{bedpe_out_fh});

  unless ($count || ($chunk > 1)) {
	  warn "NO ENTRIES PRINTED.\n";
  }

  return($count);
}
#-----------------------------------------------------------------------------#

=head2 parse_header

  Arg (1)    : brassI file_handle
  Example    : $object->parse_header();
  Description: parses brassI file header line for a list of sample names
  Return     : reference to an array of samples

=cut


# go through '#' header lines
sub _parse_header {
    my $self = shift;
    my $fh = shift;

    my $header = '';
    my $sample_count = 0;
    my $samples_found = 0;
    my @samples = ();
    while (my $line = <$fh>) {
	next unless ($line =~ /\S/);

	if ($line =~ /^#/) { $header .= $line; }
	else               { last; } # end of header reached

	if ($line =~ /\#\s*NSAMPLES\s+(\d+)/) { $sample_count = $1; }
	elsif ($line =~ /\#\s*SAMPLE\s+\d+\s+(\S+)/) {
	    push @samples, $1;
	    $samples_found++;
	    last if ($samples_found && ($sample_count == $samples_found));
	}
    }
    $self->samples(\@samples);
    $self->header($header);
    return(\@samples, $header);
}
#-----------------------------------------------------------------------------#
sub _add_row {
    my ($self, $line, $rearrangement_id) = @_;
# 1 +       826169     826234     1 -       825741     825801    5       0       IL12_4199:8:68:11113:18351;IL12_4199:5:15:6783:1698;IL12_4199:2:108:3821:17066;IL12_4199:2:47:3852:4706;IL12_4199:3:102:19210:9065      -       REPEAT  AluS, LTR41

    my ($passed_checks, $repeats, $sample_read_data, $np_sample_count, $np_count, $normal_count, $tumour_count, $distance) = $self->_check_row(@$line);
    return(0) unless ($passed_checks);

    my ($chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3, @the_rest) = @$line;

    $self->_print_row($chrL,$strandL,$L5,$L3,$chrH,$strandH,$H5,$H3,$repeats,$np_sample_count,$tumour_count,$normal_count,$np_count,$distance,$sample_read_data,$rearrangement_id);
    return(1);
}

#-----------------------------------------------------------------------------#

sub _check_row {
  my ($self, $chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3, @the_rest) = @_;
  # 1 +       826169     826234     1 -       825741     825801    5       0       IL12_4199:8:68:11113:18351;IL12_4199:5:15:6783:1698;IL12_4199:2:108:3821:17066;IL12_4199:2:47:3852:4706;IL12_4199:3:102:19210:9065      -       REPEAT  AluS, LTR41


  # skip if either end maps to decoy sequences
  my $skip_chr = 'hs37d5cs';
   return(0) if ($chrL eq $skip_chr || $chrH eq $skip_chr);

  # skip it if low and high ranges are exactly the same
  if (($chrL eq $chrH) && ($L5 == $H5) && ($L3 == $H3)) {
    print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. L and H coords the same. skip\n" if ($self->{debug});
    return(0);
  }

  # process the names lists and the repeats list
  # @the_rest contains counts for each sample, read_name lists for each sample, repeats
  my ($sample_read_data, $repeats) = $self->_process_reads_and_repeats(@the_rest);

  # skip if there are repeats and 'discard_if_repeats' has been set
  if ($self->{discard_if_repeats} && $repeats) {
    print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. Repeats. skip\n" if ($self->{debug});
    return(0);
  }

  # get the read_counts and check that at least one tumour or metastatic sample exceeds the min_tumour_count (readcount) threshold
  my ($load_check, $np_sample_count, $np_count, $normal_count, $tumour_count) = $self->_get_read_counts($sample_read_data);

  # check counts meet filtering parameters
  # Susie requested this to handle metastases too - 27/4/12
  # at least one tumour or metastatic sample exceeds the min_tumour_count (readcount) threshold
  unless ($load_check) {
    print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. load_check failed. skip\n" if ($self->{debug});
    return(0);
  }
  # Susie suggested this to limit loading too much rubbish - 25/7/13
  if ($normal_count > $self->{max_normal_count}) {
    print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. max_normal_count failed. skip\n" if ($self->{debug});
    return(0);
  }

  my $is_small_foldback = 0;
  $is_small_foldback = 1 if( $chrL eq $chrH
                          && $strandL eq $strandH
                          && (($H5 + $H3)/2) - (($L5 + $L3)/2) <= 5000); # gives the mid point of the group ranges

  my $max_np_count = $self->{max_np_count};

  # Susie suggested this to limit loading too much rubbish - 25/7/13
  # the +2 for small foldback by Yilong
  $max_np_count+=2 if($is_small_foldback);

  if ($np_count > $max_np_count) {
    if ($self->{debug}) {
      if($is_small_foldback) {
        print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. max_np_count+2 failed (small foldback). skip\n";
      }
      else {
        print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. max_np_count failed. skip\n";
      }
    }
    return(0);
  }
  # Susie suggested this to limit loading too much rubbish - 27/4/12
  if ($np_sample_count > $self->{max_np_sample_count}) {
    print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. max_np_sample_count failed. skip\n" if ($self->{debug});
    return(0);
  }

  # Do basic sanity checks:
  # check that second coordinates is lower than first
  if (($chrL eq $chrH) && ($L5 > $H5)) {
    warn "\nERROR: ".join("\t",$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3)."\n";
    die "ERROR: L higher than H coord\n";
  }

  # check distance threshold exceeded
  my $distance = undef;
  if ($chrL eq $chrH) {
    $distance = $H5 - $L5;
    unless ($distance > $self->{distance_threshold}) {
      print "$chrL, $strandL, $L5, $L3, $chrH, $strandH, $H5, $H3. distance_check failed. skip\n" if ($self->{debug});
      return(0);
    }
  }

  my $passed_checks = 1;

  return($passed_checks, $repeats, $sample_read_data, $np_sample_count, $np_count, $normal_count, $tumour_count, $distance);
}
#-----------------------------------------------------------------------------#
sub _process_reads_and_repeats {
  my $self = shift;
  my @the_rest = @_;

  my $sample_read_data = {};

  # process the names list counts
  my $sample_list = $self->samples();

  foreach my $sample(@$sample_list) {
    my $sample_type = $self->sample_type($sample);
    $sample_read_data->{$sample}->{count} = shift @the_rest;
  }
  # process the repeats
  my $repeats = shift(@the_rest);
  if ($repeats && !($repeats =~ /\./)) {
    $repeats =~ s/\s+$//; # make sure there are no spaces at the end of the list
    $repeats =~ s/\s+/, /g; # make sure each repeat is in a comma separated list
#    $repeats =~ s/[^,]\s+/, /; # make sure each repeat is in a comma separated list
#    $repeats =~ s/,(\S+)/, $1/; # ...and that each comma is followed by a space
  }
  else { $repeats = ''; }

  # process the names lists
  foreach my $sample2(@$sample_list) {
    my $read_list = shift @the_rest;

    next if($read_list =~ m/^[\.\-]$/);

    # make sure the same readname does not appear twice in the readlist
    my %reads_hash;
    if($read_list) {
      for(split /;/, $read_list) {
        $reads_hash{$_} = 1;
      }
    }
    my @read_list = (sort {$a cmp $b} keys %reads_hash);
    $sample_read_data->{$sample2}->{read_names} = \@read_list;

  }
  return($sample_read_data, $repeats);
}
#-----------------------------------------------------------------------------#
    # work out the tumour/normal/metastasis counts
sub _get_read_counts {
  my ($self, $sample_read_data) = @_;

  my $np_sample_count = 0;
  my $np_count = 0;
  my $normal_count = 0;
  my $tumour_count = 0;
  my $metastasis_count = 0; # may be more than one mestastatic sample

  my $min_tumour_count = $self->{min_tumour_count};
  my $load_check = 0; # check that we have a tumour or metastatic sample which has readcount > min_tumour_count

  my $sample_list = $self->samples();

  foreach my $sample(@$sample_list) {
    next unless (defined $sample_read_data->{$sample}->{read_names});
    my $read_count = scalar(@{$sample_read_data->{$sample}->{read_names}});
    next if($read_count == 0);

    my $sample_type = $self->sample_type($sample);

    if ($sample_type eq 'T')    {
      $tumour_count =  $read_count;
      $load_check = 1 if ($tumour_count >= $min_tumour_count);
    }
    elsif ($sample_type eq 'N') {
      $normal_count =  $read_count;
    }
    elsif ($sample_type eq 'M') {
      my $m_count = $read_count;
      $load_check = 1 if ($m_count >= $min_tumour_count);
      $metastasis_count =  $metastasis_count + $m_count;
    }
    elsif ($sample_type eq 'P') {
      $np_count =  $np_count + $read_count; # number or reads in all the normal panel samples
      $np_sample_count++; # number of samples in the normal panel containing reads in this rearrangement group
    }
  }

  return($load_check, $np_sample_count, $np_count, $normal_count, $tumour_count);
}
#-----------------------------------------------------------------------------#
sub _print_row {
    my $self = shift;
    my($chrL,$strandL,$L5,$L3,$chrH,$strandH,$H5,$H3,$repeats,$np_sample_count,$tumour_count,$normal_count,$np_count,$distance,$sample_read_data,$rearrangement_id) = @_;

    my $bedpe_fh = $self->{bedpe_out_fh};
    my $sample_list = $self->samples();

    unless ($distance) {
	if ($chrH eq $chrL) { $distance = 0; }
	else { $distance = -1; }
    }

    # values required for bedpe format (zero referenced start positions)
    my $start1 = $L5 - 1;
    my $start2 = $H5 - 1;
    my $total_read_count = $tumour_count + $normal_count + $np_count; # to put in the 'score' field of bedpe format

    # go through all the samples in the list, print the list of read_name and other data for each
    foreach my $sample(@$sample_list) {

	my $sample_type = $self->sample_type($sample);

	next unless ($sample_type eq 'T');

	my $count = 0;
	my $names = '';
	if ($sample_read_data->{$sample}->{read_names} && @{$sample_read_data->{$sample}->{read_names}}) {
	    $count = scalar(@{$sample_read_data->{$sample}->{read_names}});
	    $names = join ',', @{$sample_read_data->{$sample}->{read_names}};
	}

        # output just the tumour entries in bedpe format (zero referenced start positions)
	print $bedpe_fh "$chrL\t$start1\t$L3\t$chrH\t$start2\t$H3\t$rearrangement_id\t$total_read_count\t$strandL\t$strandH\t$repeats\t$np_sample_count\t$tumour_count\t$normal_count\t$np_count\t$distance\t$sample\t$sample_type\t$names\t$count\t0\t0\t0\t0\t0\t0\n";
    }
    $self->{bedpe_out_fh} = $bedpe_fh;
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

1;
