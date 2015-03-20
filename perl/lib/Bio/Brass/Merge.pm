package Bio::Brass::Merge;

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


BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)};
};

use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use List::Util qw(any);

use Bio::Brass qw($VERSION);
use Bio::Brass::Group;

use Tabix;

use Const::Fast qw(const);

const my @EXPECTED => qw(comment_char normal_groups analysis_groups tumour);

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
  if((ref $args[0]) eq 'HASH') {
    %options = %{$args[0]};
  }
  else {
    %options = @args;
  }
  for my $key(keys %options) {
    next unless(any {$key eq $_} @EXPECTED);
    $self->{"_$key"} = $options{$key};
  }
  $self->{'normal_groups'} = Bio::Brass::Group->new(header => $self->tabix_header);
  $self->{'analysis_groups'} = Bio::Brass::Group->new(header => $self->groups_header, sample_chk => $self->{'_tumour'});

  # overide the file for the normal groups
  $self->{'normal_groups'}->update_input($self->{'_normal_groups'});
}

sub merge_records {
  my ($self, $ofh) = @_;
  my $comment = $self->{'_comment_char'};
  my $a_grp = $self->{'analysis_groups'};
  my $n_grp = $self->{'normal_groups'};
  my $brass_np = new Tabix(-data => $self->{'_normal_groups'});
  my $fh = $self->{'_analysis_fh'};
  while(my $line = <$fh>) {
    next if($line =~ m/^$comment/);

    # always clear the existing groups
    $a_grp->clear_group;
    $n_grp->clear_group;

    chomp $line;
    next unless($a_grp->new_group($line));
    # tabix always works half-open, groups file is 1 based
    my $res = $brass_np->query($a_grp->low_chr, $a_grp->low_5p - 1, $a_grp->low_3p);
    my %overlaps;
    if(defined $res->get) {
      while(my $record = $brass_np->read($res)){
        $n_grp->new_group($record);
        next if($a_grp->low_strand ne $n_grp->low_strand);
        next if($a_grp->high_strand ne $n_grp->high_strand);
        next if($a_grp->high_chr ne $n_grp->high_chr);
        next unless($a_grp->high_3p >= $n_grp->high_5p && $a_grp->high_5p <= $n_grp->high_3p);
        $overlaps{"@{$n_grp->{_loc_data}}"} = $record;
      }
    }

    my $overlap = filter_overlaps(\%overlaps, $n_grp->sample_count);
    if(defined $overlap) {
      $n_grp->new_group($overlap);
    }
    else {
      $n_grp->clear_group;
    }
    print $ofh $self->merge_record, "\n";
  }
  close $fh;
  delete $self->{'_analysis_fh'};
}

sub filter_overlaps {
  my ($overlaps, $samples) = @_;
#use Data::Dumper;
#print Dumper(\@_);
#exit;
  my $ret;
  # must return a single value
  # Choose the overlapping group with the highest support (most samples first, then reads)
  my @keys = keys %{$overlaps};
  my $count = (scalar @keys);
  return $ret if($count == 0);
  return $overlaps->{$keys[0]} if($count == 1);
  my $max_samples = 0;
  my $max_reads = 0;
  my $best_key;
  my $count_end_pos = 7 + $samples;
  for my $key(@keys) {
    my @record = split "\t", $overlaps->{$key};
    my @read_counts = (@record)[8..$count_end_pos];
    my $have_counts = 0;
    my $count_reads = 0;
    for(@read_counts) {
      if($_ > 0) {
        $have_counts++;
        $count_reads+= $_;
      }
    }
    if($have_counts >= $max_samples) {
      if($have_counts == $max_samples && $count_reads > $max_reads) {
        $max_reads = $count_reads;
        $best_key = $key;
      }
      else {
        $max_samples = $have_counts;
        $best_key = $key;
      }
    }
  }
  return $overlaps->{$best_key};
}

sub merge_record {
  my $self = shift;
  # first element is the analysis, second the normals
  # if n_grp has no group then just pad the record correctly
  my $a_grp = $self->{'analysis_groups'};
  my $n_grp = $self->{'normal_groups'};

  my (@n_read_counts, @n_read_lists);
  if($n_grp->read_counts) {
    @n_read_counts = $n_grp->read_counts;
    @n_read_lists = $n_grp->read_lists;
  }
  else {
    unless(exists $self->{'_pad_read_counts'}) {
      $self->{'_pad_read_counts'} = [(0) x $n_grp->sample_count];
      $self->{'_pad_read_lists'} = [(q{.}) x $n_grp->sample_count];
    }
    @n_read_counts = @{$self->{'_pad_read_counts'}};
    @n_read_lists = @{$self->{'_pad_read_lists'}};
  }

  return join("\t", $a_grp->loc_data,
                    $a_grp->read_counts, @n_read_counts,
                    $a_grp->repeats,
                    $a_grp->read_lists, @n_read_lists);
}

sub merge_headers {
  # analysis first, normals second
  my $self = shift;
  my $new_header = join "\n", $self->{'analysis_groups'}->info;
  $new_header .= "\n#INPUT\t";
  $new_header .= join "\n#INPUT\t", ($self->{'analysis_groups'}->inputs, $self->{'normal_groups'}->inputs);
  $new_header .= "\n#NSAMPLES\t".($self->{'analysis_groups'}->sample_count + $self->{'normal_groups'}->sample_count)."\n";
  my $sample_counter = 1;
  for my $sample($self->{'analysis_groups'}->samples, $self->{'normal_groups'}->samples) {
    $new_header .= sprintf "#SAMPLE\t%s\t%s\n", $sample_counter++, $sample;
  }
  return $new_header;
}

=item tabix_header

Reads any header attached to the Tabix indexed normal_groups file.

Requires 'comment_char' to be defined.

=cut

sub tabix_header {
  my $self = shift;
  # need a little code to read the header from the bgzipped file
  my $groups_header = q{};
  my $comment = $self->{'_comment_char'};
  my $z = IO::Uncompress::Gunzip->new($self->{'_normal_groups'}) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  while(my $line = $z->getline()) {
    last unless($line =~ m/^$comment/);
    $groups_header .= $line;
  }
  $z->close;
  return $groups_header
}

=item groups_header

Reads any header attached to the analysis_groups file.

Requires 'comment_char' to be defined.

=cut

sub groups_header {
  my $self = shift;
  # need a little code to read the header from the bgzipped file
  my $groups_header = q{};
  my $comment = $self->{'_comment_char'};
  my $fh;
  if($self->{'_analysis_groups'} eq q{-}) {
    $fh = *STDIN;
  }
  else {
    open $fh, '<', $self->{'_analysis_groups'};
  }
  while(my $line = <$fh>) {
    last unless($line =~ m/^$comment/);
    $groups_header .= $line;
  }
  $self->{'_analysis_fh'} = $fh;
  return $groups_header;
}

sub DESTROY {
  my $self = shift;
  close $self->{'_analysis_fh'} if(exists $self->{'_analysis_fh'});
}

1;
