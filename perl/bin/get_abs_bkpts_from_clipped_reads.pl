#!/usr/bin/perl

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

use warnings FATAL => 'all';
use Bio::DB::HTS;
use Getopt::Long;
use List::Util qw(first min max);
use Statistics::Basic qw(median mean);
use Const::Fast qw(const);
use Storable qw(dclone);

const my $SLOP => 3;

my $SLOP_LEFT_FOR_GETTING_READS = 900;
my $SLOP_RIGHT_FOR_GETTING_READS = 100;
my $CLIPPED_READS_NEEDED = 2;
my $fa_file = '';
my $out_file = '';
GetOptions(
  'slop_left_for_getting_reads=i' => \$SLOP_LEFT_FOR_GETTING_READS,
  'slop_right_for_getting_reads=i' => \$SLOP_RIGHT_FOR_GETTING_READS,
  'clipped_reads_needed=i' => \$CLIPPED_READS_NEEDED,
  'fasta=s' => \$fa_file,
  'out=s' => \$out_file,
);

my $bam_file = shift;
my $input = shift;

unless(-e $bam_file && -e $input) {
  die "\nUSAGE: get_abs_bkpts_from_clipped_reads.pl [options] full.bam merge_double_rgs.out\n\n";
}

my $bam = Bio::DB::HTS->new(-bam => $bam_file);
my $fai = Bio::DB::HTS::Fai->load($fa_file);

my @lines;
open my $IN, '<', $input || die $!;
while(my $line = <$IN>) {
  next if($line =~ m/^#/);
  chomp $line;
  my @F = (split /\t/, $line)[0..9];
  push @lines, \@F;
}
close $IN;

my (%plus_rg_bkpts_of_chr, %minus_rg_bkpts_of_chr);

for (@lines) {
  my @F = @{$_};
  if ($F[8] eq '+') {
    push @{$plus_rg_bkpts_of_chr{$F[0]}}, int(($F[1]+1+$F[2])/2);
  }
  else {
    push @{$minus_rg_bkpts_of_chr{$F[0]}}, int(($F[1]+1+$F[2])/2);
  }
  if ($F[9] eq '+') {
    push @{$plus_rg_bkpts_of_chr{$F[3]}}, int(($F[4]+1+$F[5])/2);
  }
  else {
    push @{$minus_rg_bkpts_of_chr{$F[3]}}, int(($F[4]+1+$F[5])/2);
  }
}

my %no_end_reads;

print STDERR "Collecting supporting read pairs for each region...\n";
my (%low_end_reads_of_rg, %high_end_reads_of_rg);
for (@lines) {
  # saves all the re-splitting reusing the originally parsed data
  my @F = @{$_};
  my $l = int(($F[1]+$F[2]+1)/2);
  my $h = int(($F[4]+$F[5]+1)/2);

  # Low end bkpt
  my @reads = collect_reads_by_region(
      $bam,
      $F[0], $l, $F[8],
      $F[3], $h, $F[9],
    );
  if (@reads == 0) {
    $no_end_reads{$F[6]}{'low'} = 1;
    warn "No low_end reads for record $F[6]\n";
  }
  $low_end_reads_of_rg{$F[6]} = [@reads];

  # High end bkpt
  @reads = collect_reads_by_region(
    $bam,
    $F[3], $h, $F[9],
    $F[0], $l, $F[8],
  );
  if (@reads == 0) {
    $no_end_reads{$F[6]}{'high'} = 1;
    warn "No high_end reads for record $F[6]\n";
  }
  $high_end_reads_of_rg{$F[6]} = [@reads];
}

print STDERR "Resolving overlapping rg end ranges...\n";


for my $i (0..($#lines - 1)) {
  for my $j (($i+1)..$#lines) {
    # saves all the re-splitting reusing the originally parsed data
    # there's some modification of values here though so clone it
    my @F1 = @{ dclone($lines[$i]) };
    my @F2 = @{ dclone($lines[$j]) };
    if ($F1[0] ne $F2[0] || $F1[3] ne $F2[3]) {
      next;
    }
    $F1[1]++;
    $F1[4]++;
    $F2[1]++;
    $F2[4]++;
    my $id1 = $F1[6];
    my $id2 = $F2[6];

    # Both ends must overlap. Otherwise just go next
    if (
      $F1[0] ne $F2[0] ||
      $F1[8] ne $F2[8] ||
      lowest_pos_of_reads($low_end_reads_of_rg{$id1}) > highest_pos_of_reads($low_end_reads_of_rg{$id2}) ||
      lowest_pos_of_reads($low_end_reads_of_rg{$id2}) > highest_pos_of_reads($low_end_reads_of_rg{$id1}) ||
      $F1[3] ne $F2[3] ||
      $F1[9] ne $F2[9] ||
      lowest_pos_of_reads($high_end_reads_of_rg{$id1}) > highest_pos_of_reads($high_end_reads_of_rg{$id2}) ||
      lowest_pos_of_reads($high_end_reads_of_rg{$id2}) > highest_pos_of_reads($high_end_reads_of_rg{$id1})
    ) {
      next;
    }

    # Now try to find out whether we should use the low or high end as 'anchoring region'
    if (abs($F1[1]/2 + $F1[2]/2 - $F2[1]/2 - $F2[2]/2) >= abs($F1[4]/2 + $F1[5]/2 - $F2[4]/2 - $F2[5]/2)) {
      # Use low end as anchor
      if ($F1[8] eq '+') {
        # The low end is '+' orientation
        if ($F1[1]/2 + $F1[2]/2 <= $F2[1]/2 + $F2[2]/2) {
          # @F1 corresponds to lower group
          my ($clip_pos, $clip_count) = mode_of_high_clip_pos_of_reads(grep { $_->end >= $F1[1]/2 + $F1[2]/2 - 300 } @{$low_end_reads_of_rg{$id1}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = highest_pos_of_reads($low_end_reads_of_rg{$id1});
          }
          $clip_pos = $clip_pos + $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->end > $clip_pos } @{$low_end_reads_of_rg{$id2}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->end } (@{$low_end_reads_of_rg{$id1}}, @{$low_end_reads_of_rg{$id2}}))->query();
            @{$low_end_reads_of_rg{$id1}} = grep { $_->end <= $median_pos } @{$low_end_reads_of_rg{$id1}};
            @{$low_end_reads_of_rg{$id2}} = grep { $_->end >  $median_pos } @{$low_end_reads_of_rg{$id2}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
              @lower_group_reads  = map { $_->query->name } grep( { $_->end <= $clip_pos } @{$low_end_reads_of_rg{$id1}});
              @higher_group_reads = map { $_->query->name } grep( { $_->end >  $clip_pos } @{$low_end_reads_of_rg{$id2}});
              @{$low_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @lower_group_reads) } @{$low_end_reads_of_rg{$id1}};
              @{$low_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @higher_group_reads) } @{$low_end_reads_of_rg{$id2}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F1[0], ($F1[1]+$F1[2])/2, $F1[8], $F1[3], ($F1[4]+$F1[5])/2, $F1[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end <= $clip_pos } @{$high_end_reads_of_rg{$id1}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F2[0], ($F2[1]+$F2[2])/2, $F2[8], $F2[3], ($F2[4]+$F2[5])/2, $F2[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end > $clip_pos } @{$high_end_reads_of_rg{$id2}};
        }
        else {
          # @F2 corresponds to lower group
          my ($clip_pos, $clip_count) = mode_of_high_clip_pos_of_reads(grep { $_->end >= $F2[1]/2 + $F2[2]/2 - 300 } @{$low_end_reads_of_rg{$id2}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = highest_pos_of_reads($low_end_reads_of_rg{$id2});
          }
          $clip_pos = $clip_pos + $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->end > $clip_pos } @{$low_end_reads_of_rg{$id1}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->end } (@{$low_end_reads_of_rg{$id1}}, @{$low_end_reads_of_rg{$id2}}))->query();
            @{$low_end_reads_of_rg{$id2}} = grep { $_->end <= $median_pos } @{$low_end_reads_of_rg{$id2}};
            @{$low_end_reads_of_rg{$id1}} = grep { $_->end >  $median_pos } @{$low_end_reads_of_rg{$id1}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
            @lower_group_reads  = map { $_->query->name } grep( { $_->end <= $clip_pos } @{$low_end_reads_of_rg{$id2}});
            @higher_group_reads = map { $_->query->name } grep( { $_->end >  $clip_pos } @{$low_end_reads_of_rg{$id1}});
            @{$low_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @lower_group_reads) } @{$low_end_reads_of_rg{$id2}};
            @{$low_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @higher_group_reads) } @{$low_end_reads_of_rg{$id1}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F2[0], ($F2[1]+$F2[2])/2, $F2[8], $F2[3], ($F2[4]+$F2[5])/2, $F2[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end <= $clip_pos } @{$high_end_reads_of_rg{$id2}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F1[0], ($F1[1]+$F1[2])/2, $F1[8], $F1[3], ($F1[4]+$F1[5])/2, $F1[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end > $clip_pos } @{$high_end_reads_of_rg{$id1}};
        }
      }
      else {
        # The low end is '-' orientation
        if ($F1[1]/2 + $F1[2]/2 >= $F2[1]/2 + $F2[2]/2) {
          # @F1 corresponds to higher group
          my ($clip_pos, $clip_count) = mode_of_low_clip_pos_of_reads(grep { $_->start <= $F1[1]/2 + $F1[2]/2 + 300 } @{$low_end_reads_of_rg{$id1}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = lowest_pos_of_reads($low_end_reads_of_rg{$id1});
          }
          $clip_pos = $clip_pos - $SLOP;

          my @higher_group_reads;
          my @lower_group_reads;
          if (!(grep { $_->start < $clip_pos } @{$low_end_reads_of_rg{$id2}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->start } (@{$low_end_reads_of_rg{$id1}}, @{$low_end_reads_of_rg{$id2}}))->query();
            @{$low_end_reads_of_rg{$id1}} = grep { $_->start >= $median_pos } @{$low_end_reads_of_rg{$id1}};
            @{$low_end_reads_of_rg{$id2}} = grep { $_->start <  $median_pos } @{$low_end_reads_of_rg{$id2}};
            $clip_pos = $median_pos if $median_pos;
          }

          @higher_group_reads  = map { $_->query->name } grep( { $_->start >= $clip_pos } @{$low_end_reads_of_rg{$id1}});
          @lower_group_reads   = map { $_->query->name } grep( { $_->start <  $clip_pos } @{$low_end_reads_of_rg{$id2}});
          @{$low_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @higher_group_reads) } @{$low_end_reads_of_rg{$id1}};
          @{$low_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @lower_group_reads) } @{$low_end_reads_of_rg{$id2}};

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F1[0], ($F1[1]+$F1[2])/2, $F1[8], $F1[3], ($F1[4]+$F1[5])/2, $F1[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start >= $clip_pos } @{$high_end_reads_of_rg{$id1}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F2[0], ($F2[1]+$F2[2])/2, $F2[8], $F2[3], ($F2[4]+$F2[5])/2, $F2[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start < $clip_pos } @{$high_end_reads_of_rg{$id2}};
        }
        else {
          # @F2 corresponds to higher group
          my ($clip_pos, $clip_count) = mode_of_low_clip_pos_of_reads(grep { $_->start <= $F2[1]/2 + $F2[2]/2 + 300 } @{$low_end_reads_of_rg{$id2}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = lowest_pos_of_reads($low_end_reads_of_rg{$id2});
          }
          $clip_pos = $clip_pos - $SLOP;

          my @higher_group_reads;
          my @lower_group_reads;
          if (!(grep { $_->start < $clip_pos } @{$low_end_reads_of_rg{$id1}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->start } (@{$low_end_reads_of_rg{$id1}}, @{$low_end_reads_of_rg{$id2}}))->query();
            @{$low_end_reads_of_rg{$id2}} = grep { $_->start >= $median_pos } @{$low_end_reads_of_rg{$id2}};
            @{$low_end_reads_of_rg{$id1}} = grep { $_->start <  $median_pos } @{$low_end_reads_of_rg{$id1}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
            @higher_group_reads  = map { $_->query->name } grep( { $_->start >= $clip_pos } @{$low_end_reads_of_rg{$id2}});
            @lower_group_reads   = map { $_->query->name } grep( { $_->start <  $clip_pos } @{$low_end_reads_of_rg{$id1}});
            @{$low_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @higher_group_reads) } @{$low_end_reads_of_rg{$id2}};
            @{$low_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @lower_group_reads) } @{$low_end_reads_of_rg{$id1}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F2[0], ($F2[1]+$F2[2])/2, $F2[8], $F2[3], ($F2[4]+$F2[5])/2, $F2[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start >= $clip_pos } @{$high_end_reads_of_rg{$id2}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F1[0], ($F1[1]+$F1[2])/2, $F1[8], $F1[3], ($F1[4]+$F1[5])/2, $F1[9], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$high_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start < $clip_pos } @{$high_end_reads_of_rg{$id1}};
        }
      }
    }
    else {
      # Use high end as anchor
      if ($F1[9] eq '+') {
        # The high end is '+' orientation
        if ($F1[4]/2 + $F1[5]/2 <= $F2[4]/2 + $F2[5]/2) {
          # @F1 corresponds to lower group
          my ($clip_pos, $clip_count) = mode_of_high_clip_pos_of_reads(grep { $_->end >= $F1[4]/2 + $F1[5]/2 - 300 } @{$high_end_reads_of_rg{$id1}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = highest_pos_of_reads($high_end_reads_of_rg{$id1});
          }
          $clip_pos = $clip_pos + $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->end > $clip_pos } @{$high_end_reads_of_rg{$id2}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->end } (@{$high_end_reads_of_rg{$id1}}, @{$high_end_reads_of_rg{$id2}}))->query();
            @{$high_end_reads_of_rg{$id1}} = grep { $_->end <= $median_pos } @{$high_end_reads_of_rg{$id1}};
            @{$high_end_reads_of_rg{$id2}} = grep { $_->end >  $median_pos } @{$high_end_reads_of_rg{$id2}};
            $clip_pos = $median_posif $median_pos;
          }
          else {
            @lower_group_reads  = map { $_->query->name } grep( { $_->end <= $clip_pos } @{$high_end_reads_of_rg{$id1}});
            @higher_group_reads = map { $_->query->name } grep( { $_->end >  $clip_pos } @{$high_end_reads_of_rg{$id2}});
            @{$high_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @lower_group_reads) } @{$high_end_reads_of_rg{$id1}};
            @{$high_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @higher_group_reads) } @{$high_end_reads_of_rg{$id2}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F1[3], ($F1[4]+$F1[5])/2, $F1[9], $F1[0], ($F1[1]+$F1[2])/2, $F1[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end <= $clip_pos } @{$low_end_reads_of_rg{$id1}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F2[3], ($F2[4]+$F2[5])/2, $F2[9], $F2[0], ($F2[1]+$F2[2])/2, $F2[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end > $clip_pos } @{$low_end_reads_of_rg{$id2}};
        }
        else {
          # @F2 corresponds to lower group
          my ($clip_pos, $clip_count) = mode_of_high_clip_pos_of_reads(grep { $_->end >= $F2[4]/2 + $F2[5]/2 - 300 } @{$high_end_reads_of_rg{$id2}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = highest_pos_of_reads($high_end_reads_of_rg{$id2});
          }
          $clip_pos = $clip_pos + $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->end > $clip_pos } @{$high_end_reads_of_rg{$id1}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->end } (@{$high_end_reads_of_rg{$id1}}, @{$high_end_reads_of_rg{$id2}}))->query();
            @{$high_end_reads_of_rg{$id2}} = grep { $_->end <= $median_pos } @{$high_end_reads_of_rg{$id2}};
            @{$high_end_reads_of_rg{$id1}} = grep { $_->end >  $median_pos } @{$high_end_reads_of_rg{$id1}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
            @lower_group_reads  = map { $_->query->name } grep( { $_->end <= $clip_pos } @{$high_end_reads_of_rg{$id2}});
            @higher_group_reads = map { $_->query->name } grep( { $_->end >  $clip_pos } @{$high_end_reads_of_rg{$id1}});
            @{$high_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @lower_group_reads) } @{$high_end_reads_of_rg{$id2}};
            @{$high_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @higher_group_reads) } @{$high_end_reads_of_rg{$id1}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F2[3], ($F2[4]+$F2[5])/2, $F2[9], $F2[0], ($F2[1]+$F2[2])/2, $F2[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end <= $clip_pos } @{$low_end_reads_of_rg{$id2}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F1[3], ($F1[4]+$F1[5])/2, $F1[9], $F1[0], ($F1[1]+$F1[2])/2, $F1[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->end > $clip_pos } @{$low_end_reads_of_rg{$id1}};
        }
      }
      else {
        # The high end is '-' orientation
        if ($F1[4]/2 + $F1[5]/2 >= $F2[4]/2 + $F2[5]/2) {
          # @F1 corresponds to higher group
          my ($clip_pos, $clip_count) = mode_of_low_clip_pos_of_reads(grep { $_->start <= $F1[4]/2 + $F1[5]/2 + 300 } @{$high_end_reads_of_rg{$id1}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = lowest_pos_of_reads($high_end_reads_of_rg{$id1});
          }
          $clip_pos = $clip_pos - $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->start < $clip_pos } @{$high_end_reads_of_rg{$id2}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->start } (@{$high_end_reads_of_rg{$id1}}, @{$high_end_reads_of_rg{$id2}}))->query();
            @{$high_end_reads_of_rg{$id1}} = grep { $_->start >= $median_pos } @{$high_end_reads_of_rg{$id1}};
            @{$high_end_reads_of_rg{$id2}} = grep { $_->start <  $median_pos } @{$high_end_reads_of_rg{$id2}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
            @higher_group_reads = map { $_->query->name } grep( { $_->start >=  $clip_pos } @{$high_end_reads_of_rg{$id1}});
            @lower_group_reads  = map { $_->query->name } grep( { $_->start < $clip_pos } @{$high_end_reads_of_rg{$id2}});
            @{$high_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @higher_group_reads) } @{$high_end_reads_of_rg{$id1}};
            @{$high_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @lower_group_reads) } @{$high_end_reads_of_rg{$id2}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F1[3], ($F1[4]+$F1[5])/2, $F1[9], $F1[0], ($F1[1]+$F1[2])/2, $F1[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start >= $clip_pos } @{$low_end_reads_of_rg{$id1}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F2[3], ($F2[4]+$F2[5])/2, $F2[9], $F2[0], ($F2[1]+$F2[2])/2, $F2[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start < $clip_pos } @{$low_end_reads_of_rg{$id2}};
        }
        else {
          # @F2 corresponds to higher group
          my ($clip_pos, $clip_count) = mode_of_low_clip_pos_of_reads(grep { $_->start <= $F2[4]/2 + $F2[5]/2 + 300 } @{$high_end_reads_of_rg{$id2}});
          if ($clip_count < $CLIPPED_READS_NEEDED) {
            $clip_pos = lowest_pos_of_reads($high_end_reads_of_rg{$id2});
          }
          $clip_pos = $clip_pos - $SLOP;

          my @lower_group_reads;
          my @higher_group_reads;
          if (!(grep { $_->start < $clip_pos } @{$high_end_reads_of_rg{$id1}})) {
            # die;
            warn("Could not properly separate the reads between the low ends of rearrangements $id1 and $id2!\n");
            my $median_pos = median(map { $_->start } (@{$high_end_reads_of_rg{$id1}}, @{$high_end_reads_of_rg{$id2}}))->query();
            @{$high_end_reads_of_rg{$id2}} = grep { $_->start >= $median_pos } @{$high_end_reads_of_rg{$id2}};
            @{$high_end_reads_of_rg{$id1}} = grep { $_->start <  $median_pos } @{$high_end_reads_of_rg{$id1}};
            $clip_pos = $median_pos if $median_pos;
          }
          else {
            @higher_group_reads = map { $_->query->name } grep( { $_->start >=  $clip_pos } @{$high_end_reads_of_rg{$id2}});
            @lower_group_reads  = map { $_->query->name } grep( { $_->start < $clip_pos } @{$high_end_reads_of_rg{$id1}});
            @{$high_end_reads_of_rg{$id2}} = grep { within_array($_->query->name, @higher_group_reads) } @{$high_end_reads_of_rg{$id2}};
            @{$high_end_reads_of_rg{$id1}} = grep { within_array($_->query->name, @lower_group_reads) } @{$high_end_reads_of_rg{$id1}};
          }

          my %mate_of_read;
          for (collect_reads_by_region($bam, $F2[3], ($F2[4]+$F2[5])/2, $F2[9], $F2[0], ($F2[1]+$F2[2])/2, $F2[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id2}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start >= $clip_pos } @{$low_end_reads_of_rg{$id2}};

          undef %mate_of_read;
          for (collect_reads_by_region($bam, $F1[3], ($F1[4]+$F1[5])/2, $F1[9], $F1[0], ($F1[1]+$F1[2])/2, $F1[8], 1)) {
            $mate_of_read{$_->query->name} = $_;
          }
          @{$low_end_reads_of_rg{$id1}} = grep { defined($mate_of_read{$_->query->name}) && $mate_of_read{$_->query->name}->start < $clip_pos } @{$low_end_reads_of_rg{$id1}};
        }
      }
    }
  }
}

print STDERR "Outputting read position ranges...\n";

my $outfh;
if($out_file eq q{}) {
  $outfh = *STDOUT;
}
else {
  open $outfh, '>', $out_file || die $!;
}

for (@lines) {
  # saves all the re-splitting reusing the originally parsed data
  my @F = @{$_};

  my ($low_end_low_clip_pos, $low_end_low_clip_count) = mode_of_low_clip_pos_of_reads(@{$low_end_reads_of_rg{$F[6]}});
  my ($low_end_high_clip_pos, $low_end_high_clip_count) = mode_of_high_clip_pos_of_reads(@{$low_end_reads_of_rg{$F[6]}});
  my $low_end_bkpt;
  if ($F[8] eq '+') {
    $low_end_bkpt = ($low_end_high_clip_count >= $CLIPPED_READS_NEEDED ? $low_end_high_clip_pos : highest_pos_of_reads($low_end_reads_of_rg{$F[6]}));
  }
  else {
    $low_end_bkpt = ($low_end_low_clip_count >= $CLIPPED_READS_NEEDED ? $low_end_low_clip_pos : lowest_pos_of_reads($low_end_reads_of_rg{$F[6]}));
  }

  my ($high_end_low_clip_pos, $high_end_low_clip_count) = mode_of_low_clip_pos_of_reads(@{$high_end_reads_of_rg{$F[6]}});
  my ($high_end_high_clip_pos, $high_end_high_clip_count) = mode_of_high_clip_pos_of_reads(@{$high_end_reads_of_rg{$F[6]}});
  my $high_end_bkpt;
  if ($F[9] eq '+') {
    $high_end_bkpt = ($high_end_high_clip_count >= $CLIPPED_READS_NEEDED ? $high_end_high_clip_pos : highest_pos_of_reads($high_end_reads_of_rg{$F[6]}));
  }
  else {
    $high_end_bkpt = ($high_end_low_clip_count >= $CLIPPED_READS_NEEDED ? $high_end_low_clip_pos : lowest_pos_of_reads($high_end_reads_of_rg{$F[6]}));
  }

  my ($low_end_low_clip, $low_end_high_clip, $low_end_read_uniq, $low_end_read_count);
  my ($high_end_low_clip, $high_end_high_clip, $high_end_read_uniq, $high_end_read_count);

  if(exists $no_end_reads{$F[6]}{'low'}) {
    $low_end_bkpt = int abs (($F[1] + $F[2])/2);
    $low_end_low_clip = $low_end_bkpt.' (0)';
    $low_end_high_clip = $low_end_bkpt.' (0)';
    $low_end_read_uniq = q{};
    $low_end_read_count = 0;
  }
  else {
    $low_end_low_clip = ($low_end_low_clip_count   >= $CLIPPED_READS_NEEDED ? "$low_end_low_clip_pos ($low_end_low_clip_count)"   : lowest_pos_of_reads($low_end_reads_of_rg{$F[6]}) . ' (0)');
    $low_end_high_clip = ($low_end_high_clip_count  >= $CLIPPED_READS_NEEDED ? "$low_end_high_clip_pos ($low_end_high_clip_count)"   : highest_pos_of_reads($low_end_reads_of_rg{$F[6]}) . ' (0)');
    my @tmpReads = unique(map { $_->query->name } @{$low_end_reads_of_rg{$F[6]}});
    $low_end_read_uniq = join(',', @tmpReads);
    $low_end_read_count = scalar(@tmpReads);
  }

  if(exists $no_end_reads{$F[6]}{'high'}) {
    $high_end_bkpt = int abs (($F[4] + $F[5])/2);
    $high_end_low_clip = $high_end_bkpt.' (0)';
    $high_end_high_clip = $high_end_bkpt.' (0)';
    $high_end_read_uniq = q{};
    $high_end_read_count = 0;
  }
  else {
    $high_end_low_clip = ($high_end_low_clip_count  >= $CLIPPED_READS_NEEDED ? "$high_end_low_clip_pos ($high_end_low_clip_count)"   : lowest_pos_of_reads($high_end_reads_of_rg{$F[6]}) . ' (0)');
    $high_end_high_clip = ($high_end_high_clip_count >= $CLIPPED_READS_NEEDED ? "$high_end_high_clip_pos ($high_end_high_clip_count)" : highest_pos_of_reads($high_end_reads_of_rg{$F[6]}) . ' (0)');
    my @tmpReads = unique(map { $_->query->name } @{$high_end_reads_of_rg{$F[6]}});
    $high_end_read_uniq = join(',', @tmpReads);
    $high_end_read_count = scalar(@tmpReads);
  }

  print $outfh join(
    "\t",
    @F,
    $low_end_bkpt,
    $high_end_bkpt,
    $low_end_low_clip,
    $low_end_high_clip,
    $high_end_low_clip,
    $high_end_high_clip,
    $low_end_read_uniq,
    $low_end_read_count,
    ($F[8] eq '+' ? $low_end_high_clip_count : $low_end_low_clip_count),
    $high_end_read_uniq,
    $high_end_read_count,
    ($F[9] eq '+' ? $high_end_high_clip_count : $high_end_low_clip_count),
  ),"\n";
}
close $outfh if($out_file ne q{});

sub within_array {
  my $el = shift;
  return grep { $_ eq $el } @_;
}

sub collect_reads_by_region {
  my($bam, $chr, $pos, $dir, $mate_chr, $mate_pos, $mate_dir, $ignore_mate_pos) = @_;
  my ($start, $end, $mate_start, $mate_end);
  # Adjust regions of expected read positions based on BEDPE intervals
  if ($dir eq '+') {
    $start = $pos - $SLOP_LEFT_FOR_GETTING_READS;
    $end   = $pos + $SLOP_RIGHT_FOR_GETTING_READS;
  }
  else {
    $start = $pos - $SLOP_RIGHT_FOR_GETTING_READS;
    $end   = $pos + $SLOP_LEFT_FOR_GETTING_READS;
  }
  if ($mate_dir eq '+') {
    $mate_start = $mate_pos - $SLOP_LEFT_FOR_GETTING_READS;
    $mate_end   = $mate_pos + $SLOP_RIGHT_FOR_GETTING_READS;
  }
  else {
    $mate_start = $mate_pos - $SLOP_RIGHT_FOR_GETTING_READS;
    $mate_end   = $mate_pos + $SLOP_LEFT_FOR_GETTING_READS;
  }

  $start = int($start);
  $end = int($end);
  $mate_start = int($mate_start);
  $mate_end = int($mate_end);

  # Change strands to Perl-Sam format
  $dir = ($dir eq '+' ? 1 : -1);
  $mate_dir = ($mate_dir eq '+' ? 1 : -1);

  my @aln = $bam->get_features_by_location($chr, $start, $end);
  @aln = grep {
    # !is_supplementary($_->flag) &&
    $_->get_tag_values('FLAGS') =~ /PAIRED/ &&
    $_->get_tag_values('FLAGS') !~ /UNMAPPED/ &&
    $_->strand eq $dir &&
    (
      ( $ignore_mate_pos && $_->get_tag_values('FLAGS') =~ /MAP_PAIR/ ) ||
      ( $chr ne $mate_chr || $dir != $mate_dir ) ||
      ( $pos <= $mate_pos && $_->start <= $_->mate_start ) ||
      ( $pos  > $mate_pos && $_->start >= $_->mate_start )
    ) &&
    (
      $ignore_mate_pos ||
      (
        $_->mate_seq_id eq $mate_chr &&
        $_->mstrand   eq $mate_dir &&
        $_->mate_start  <= $mate_end &&
        $_->mate_start + $_->length - 1 >= $mate_start
      )
    )
  } @aln;
  return @aln;
}

sub lowest_pos_of_reads {
  my $r = shift @_;
  return min(map { $_->start } @{ $r }) if(@{$r} > 0);
  return 0;
}

sub mode_of_low_clip_pos_of_reads {
  my @start_clip_reads = grep { $_->cigar_str =~ /^\d+[HS]/ } @_;
  if (@start_clip_reads == 0) {
    return(undef, 0);
  }
  else {
    my %count_of_clip_pos = ();
    for (@start_clip_reads) {
      $count_of_clip_pos{$_->start}++;
    }
    my $clip_pos = (sort { $count_of_clip_pos{$b} <=> $count_of_clip_pos{$a} || $b <=> $a } keys(%count_of_clip_pos))[0];
    return(
      $clip_pos,
      $count_of_clip_pos{$clip_pos}
    );
  }
}

sub highest_pos_of_reads {
  my $r = shift @_;
  return max(map { $_->end } @{ $r }) if(@{$r} > 0);
  return 0;
}

sub mode_of_high_clip_pos_of_reads {
  my @end_clip_reads = grep { $_->cigar_str =~ /\d+[HS]$/ } @_;
  if (@end_clip_reads == 0) {
    return(undef, 0);
  }
  else {
    my %count_of_clip_pos = ();
    for (@end_clip_reads) {
      $count_of_clip_pos{$_->end}++;
    }
    my $clip_pos = (sort { $count_of_clip_pos{$b} <=> $count_of_clip_pos{$a} || $a <=> $b } keys(%count_of_clip_pos))[0];
    return(
      $clip_pos,
      $count_of_clip_pos{$clip_pos}
    );
  }
}

sub unique {
  my %seen;
  for (@_) {
    $seen{$_} = 1;
  }
  return sort {$a cmp $b} keys %seen;
}
