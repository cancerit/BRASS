package Bio::Brass::Group;

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


use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use List::Util qw(any);

use Bio::Brass qw($VERSION);

use Const::Fast qw(const);

=head

This is a data object
new() is created with the header text
new_group() is passed the data for a group

2 of these objects are passed to MergeGroup->new_header to create a new header
2 of these objects (with correct new_group entries) are passed to MergeGroup->new_record to merge the records into a singel element

Control of the matching is the responsibility of the wrapping code (Bio::Brass:Merge)

=cut

const my @EXPECTED => qw(header);
const my $CHR => 0;
const my $STRAND => 1;
const my $START => 2;
const my $END => 3;
const my $HIGH_OFFSET => 4;
const my $MIN_GAP => 200;

sub new {
  my ($class, @args) = @_;
  my $self = { };
  bless $self, $class;
  $self->_init(@args);
  return $self;
}

sub _init {
  my ($self, @args) = @_;
  my %options;
  if(ref $args[0] eq 'HASH') {
    %options = %{$args[0]};
  }
  else {
    %options = @args;
  }
  for my $key(keys %options) {
    next unless(any { $key eq $_} @EXPECTED);
    $self->{"_$key"} = $options{$key};
  }
  $self->parse_header($options{'sample_chk'});
}

sub parse_header {
  my ($self, $sample_chk) = @_;
  my @lines = split "\n", $self->{'_header'};
  my @inputs;
  my @samples;
  my @info;
  for my $line (@lines) {
    if($line =~ m/^#\s/) {
      push @info, $line;
      next;
    }
    if($line =~ m/^#INPUT\s(.+)/) {
      push @inputs, $1;
      next;
    }
    if($line =~ m/^#SAMPLE\s[[:digit:]]+\s+(.+)/) {
      push @samples, $1;
    }
  }

  if(defined $sample_chk) {
    # need to store the index of the tumour for later filtering
    for my $i(0..((scalar @samples) -1)) {
      if($samples[$i] eq $sample_chk) {
        $self->{'tumour_idx'} = $i;
        last;
      }
    }
    die "ERROR: Specified tumour '$sample_chk' not found in input\n"unless(exists $self->{'tumour_idx'});
  }
  if(@info) {
    splice @info, 1, 0, '#'; # put the blank line back
  }
#  else {
#    push @info, '#';
#  }
  $self->{'_inputs'} = \@inputs;
  $self->{'_samples'} = \@samples;
  $self->{'_sample_count'} = (scalar @samples);
  $self->{'_info'} = \@info;
  1;
}

sub clear_group {
  my $self = shift;
  delete $self->{'_loc_data'};
  delete $self->{'_read_counts'};
  delete $self->{'_repeats'};
  delete $self->{'_read_lists'};
  1;
}

sub new_group {
  #lChr lStr l5p l3p hChr hStr h5p h3p [sample counts] repeat(s) [sample reads]
  # process it backwards to save messy handling
  my ($self, @record) = @_;
  my $ref = ref $record[0];
  if($ref eq 'ARRAY') {
    @record = @{$record[0]};
  }
  elsif($ref eq q{}) {
    @record = split "\t", $record[0];
  }
  my $sample_count = $self->sample_count;
  my @read_lists = splice(@record, -$sample_count);

  # we can abandon construction if sample_chk is in use and the specified sample has no reads
  if(exists $self->{'tumour_idx'} && $read_lists[$self->{'tumour_idx'}] eq q{.}) {
    return 0;
  }


  my $repeats = pop @record;
  my @read_counts = splice(@record, -$sample_count);

  # ignore close events
  if(exists $self->{'tumour_idx'}) {
    my ($lChr, $lStr, $l5p, $l3p, $hChr, $hStr, $h5p, $h3p) = @record;
    if($lChr eq $hChr && $lStr eq q{+} && $hStr eq q{-} && ($h5p - $l3p) <= $MIN_GAP) {
      return 0;
    }
  }

  #my ($lChr $lStr $l5p $l3p $hChr $hStr $h5p $h3p) = @record;
  # in the order it should be reconstructed
  $self->{'_loc_data'} = \@record;
  $self->{'_read_counts'} = \@read_counts;
  $self->{'_repeats'} = $repeats;
  $self->{'_read_lists'} = \@read_lists;

  1;
}

sub low_chr {
  shift->{'_loc_data'}->[$CHR];
}

sub low_strand {
  shift->{'_loc_data'}->[$STRAND];
}

sub low_5p {
  shift->{'_loc_data'}->[$START];
}

sub low_3p {
  shift->{'_loc_data'}->[$END];
}

sub high_chr {
  shift->{'_loc_data'}->[$CHR+$HIGH_OFFSET];
}

sub high_strand {
  shift->{'_loc_data'}->[$STRAND+$HIGH_OFFSET];
}

sub high_5p {
  shift->{'_loc_data'}->[$START+$HIGH_OFFSET];
}

sub high_3p {
  shift->{'_loc_data'}->[$END+$HIGH_OFFSET];
}

sub update_input {
  my ($self, $new_input) = @_;
  # can only be applied if a single input, i.e. the header of brass_np.bam file can be replaces with brass_np.groups
  die 'update_input is only possible if input count = 1' unless(scalar @{$self->inputs} == 1);
  $self->{'_inputs'} = [$new_input];
}

sub inputs {
  my $self = shift;
  return wantarray ? @{$self->{'_inputs'}} : $self->{'_inputs'};
}

sub samples {
  my $self = shift;
  return wantarray ? @{$self->{'_samples'}} : $self->{'_samples'};
}

sub sample_count {
  shift->{'_sample_count'};
}

sub loc_data {
  my $self = shift;
  return wantarray ? @{$self->{'_loc_data'}} : $self->{'_loc_data'};
}

sub read_counts {
  my $self = shift;
  return wantarray ? @{$self->{'_read_counts'}} : $self->{'_read_counts'};
}

sub repeats {
  shift->{'_repeats'};
}

sub read_lists {
  my $self = shift;
  return wantarray ? @{$self->{'_read_lists'}} : $self->{'_read_lists'};
}

sub info {
  my $self = shift;
  return wantarray ? @{$self->{'_info'}} : $self->{'_info'};
}

1;
