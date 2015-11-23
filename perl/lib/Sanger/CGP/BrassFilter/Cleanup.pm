package Sanger::CGP::BrassFilter::Cleanup;

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
use List::Util qw(min max sum);
use Data::Dumper;

use Bio::Brass;
our $VERSION = Bio::Brass->VERSION;

use Const::Fast qw(const);

const my %DISPATCH => ('chr1' => \&die_if_diff,
                'start1'  => \&l_min,
                'end1'    => \&l_max,
                'start2'  => \&l_min,
                'chr2' => \&die_if_diff,
                'end2'    => \&l_max,
                'id/name'      => \&l_join,
                'brass_score' => \&l_read_count,
                'strand1' => \&die_if_diff,
                'strand2' => \&die_if_diff,
                'repeats' => \&merge_repeats,
                'np_sample_count' => \&l_sum,
                'tumour_count' => \&l_read_count,
                'normal_count' => \&l_sum,
                'np_count' => \&l_sum,
                'bkpt_distance' => \&correct_bp,
                'sample' => \&uniq_list,
                'sample_type' => \&uniq_list,
                'names' => \&uniq_list,
                'count' => \&sum_counts,
                'bal_trans' => \&die_if_diff,
                'inv' => \&l_join,
                'occL' => \&no_op,
                'occH' => \&no_op,
                'copynumber_flag' => \&die_if_diff,
                'range_blat' => \&l_join,);

sub new {
  my ($class, %args) = @_;
  my $self = {};
  bless $self,$class;

  $self->{'infile'} = $args{-infile} if(defined $args{-infile});
  $self->{'outfile'} = $args{-outfile} if(defined $args{-outfile});

  return $self;
}

sub process {
  my $self = shift;
  my @headers;

  my $in_file = $self->{'infile'};
  my $out_file = $self->{'outfile'};

  open my $IN, '<', $in_file;
  my $last_header;
  my (%loc_by_name, %name_by_loc);
  my @header_keys;
  my %counts;
  my %end_count;
  my @pre_clean;

  while(my $line = <$IN>) {
    if($line =~ m/^#/) {
      chomp $line;
      push @headers, $line;
      next;
    }
    unless(defined $last_header) {
      $last_header = $headers[-1];
      my @lookups = header_map($last_header);
      %loc_by_name = %{$lookups[0]};
      %name_by_loc = %{$lookups[1]};
      @header_keys = @{$lookups[2]};
    }
    chomp $line;
    $counts{'Groups Total'}++;
    my %group;
    my @elements = split /\t/, $line;
    for my $i(0..((scalar @elements)-1)) {
      my $name = $name_by_loc{$i};
      $elements[$i] =~ s/\s+//g if($name eq 'repeats');
      $group{ $name_by_loc{$i} } = $elements[$i];
    }
    $group{'_orig'} = $line;
    push @pre_clean,  \%group;
    $end_count{"$group{chr1}:$group{start1}-$group{end1}:$group{strand1}"}++;
    $end_count{"$group{chr2}:$group{start2}-$group{end2}:$group{strand2}"}++;
  }
  close $IN;

#for(sort keys %end_count){
#  warn $_.' : '.$end_count{$_}."\n" if $end_count{$_} > 10;
#}

  my $grps_before_clean = scalar @pre_clean;
  my $end_occ_count = 0;
  my @groups;
  while(my $group = shift @pre_clean) {
    next if($end_count{ "$group->{chr1}:$group->{start1}-$group->{end1}:$group->{strand1}" } > 10);
    next if($end_count{ "$group->{chr2}:$group->{start2}-$group->{end2}:$group->{strand2}" } > 10);
    push @groups, $group;
    $end_occ_count++;
  }

  $counts{'Discard - EndOcc'} = $counts{'Groups Total'} - $end_occ_count;

  my $col_max = (scalar (keys %name_by_loc))-1;
  open my $OUT, '>', $out_file;
  print $OUT join("\n", @headers), "\n";
  for my $group(@groups) {
    for my $i(0..$col_max) {
      print $OUT $group->{ $name_by_loc{$i} };
      print $OUT "\t" if($i != $col_max);
    }
    print $OUT "\n";
  }
  close $OUT;

  for(sort keys %counts) {
    warn "$_:\t$counts{$_}\n";
  }
}

sub occX {
  my $groups = shift;
  for my $this(@{$groups}) {
    my ($low, $high) = (0,0);
    for my $that(@{$groups}) {
      # I count myself too
      $low++ if(($this->{'start1'}-500) <= $that->{'end1'} && ($this->{'end1'}+500) >= $that->{'start1'});
      $high++ if(($this->{'start2'}-500) <= $that->{'end2'} && ($this->{'end2'}+500) >= $that->{'start2'});
    }
    $this->{'occL'} = $low;
    $this->{'occH'} = $high;
  }
}

sub filter_groups {
  my ($groups, $counts) = @_;
  my @new_groups;
  while(my $group = shift @{$groups}) {
    if(discard_inv($group) == 1) {
      $counts->{'Discard - Small Inv'}++;
      next;
    }
    if(discard_smallDel($group) == 1) {
      $counts->{'Discard - Small del'}++;
      next;
    }
    push @new_groups, $group;
  }
  return \@new_groups;
}

sub merge_groups {
  my ($groups, $header_keys) = @_;
  my @new_groups;
  while (my $this = shift @{$groups}) {
    # munge the coords for this event if not done in usage
    unless(exists $this->{'_start1'}) {
      $this->{'_start1'} = $this->{'start1'}-250;
      $this->{'_start2'} = $this->{'start2'}-250;
      $this->{'_end1'} = $this->{'end1'}+250;
      $this->{'_end2'} = $this->{'end2'}+250;
    }
    my $merged = 0;
    for my $that(@{$groups}) {
      # will only ever be called once per group as modifying the primary reference
      unless(exists $that->{'_start1'}) {
        $that->{'_start1'} = $that->{'start1'}-250;
        $that->{'_start2'} = $that->{'start2'}-250;
        $that->{'_end1'} = $that->{'end1'}+250;
        $that->{'_end2'} = $that->{'end2'}+250;
      }

      my $a_s1 = $this->{'_start1'};
      my $a_e1 = $this->{'_end1'};
      my $a_s2 = $this->{'_start2'};
      my $a_e2 = $this->{'_end2'};

      my $b_s1 = $that->{'_start1'};
      my $b_e1 = $that->{'_end1'};
      my $b_s2 = $that->{'_start2'};
      my $b_e2 = $that->{'_end2'};

      if(($a_s1 <= $b_e1 && $a_e1 >= $b_s1) && ($a_s2 <= $b_e2 && $a_e2 >= $b_s2)) {
        merge_event($this, $that, $header_keys);
        $merged = 1;
        # clear these as they will now need a different range
        delete $that->{'_start1'};
        delete $that->{'_start2'};
        delete $that->{'_end1'};
        delete $that->{'_end2'};
        last; # can only merge once in a loop
      }
    }

    push @new_groups, $this unless($merged);
  }
  return \@new_groups;
}

sub discard_smallDel {
  my $group = shift;
  # must be same chr, different strand to filter
  return 0 if($group->{'chr1'} ne $group->{'chr2'});
  return 0 if($group->{'strand1'} ne '+' || $group->{'strand2'} ne '-');
  my $dist = $group->{'start2'} - $group->{'end1'};
  return 1 if($dist < 500);
  return 0;
}

sub discard_inv {
  my $group = shift;
  # is it an inversion
  return 0 if(($group->{'chr1'} ne $group->{'chr2'}) || ($group->{'strand1'} ne $group->{'strand2'}));

  my $t_count = $group->{'tumour_count'};

  # some fail fast checks
  return 1 if($t_count < 4);
  return 0 if($t_count >= 10);

  my $dist = 1_000_000_000_000; # I just need a big number to start

  for my $g1 ($group->{'start1'}, $group->{'end1'}) {
    for my $g2 ($group->{'start2'}, $group->{'end2'}) {
      my $abs_dist = abs $g1 - $g2;
      $dist = $abs_dist if($abs_dist < $dist);
    }
  }

  return 0 if($dist >= 5000 && $t_count >= 4);
  return 0 if($dist < 5000 && $t_count >= 10);
  return 1;
}

sub header_map {
  my $header_line = shift;
  $header_line =~ s/^#\s+//;
  my @elements = split /\t/, $header_line;
  my (%loc_by_name, %name_by_loc);
  my @header_keys;
  my @to_end;
  my $count = 0;
  for my $e(@elements) {
    if($e =~ m/_count$/ || $e eq 'brass_score') {
      unshift @to_end, $e;
    } elsif($e eq 'count' || $e eq 'bkpt_distance') {
      # have to be done last
      push @to_end, $e;
    } else {
      push @header_keys, $e;
    }
    $loc_by_name{$e} = $count;
    $name_by_loc{$count} = $e;
    $count++;
  }
  push @header_keys, @to_end;
  return (\%loc_by_name, \%name_by_loc, \@header_keys);
}

sub uniq_list {
  my ($key, $a, $b) = @_;
  my %d = map {$_ => 1} (split /,/, $a->{$key}),(split /,/, $b->{$key});
  return (join ',', keys %d);
}

sub l_min {
  my ($key, $a, $b) = @_;
  return min($a->{$key}, $b->{$key});
}

sub l_max {
  my ($key, $a, $b) = @_;
  return max($a->{$key}, $b->{$key});
}

sub l_sum {
  my ($key, $a, $b) = @_;
  return sum($a->{$key}, $b->{$key});
}

sub l_read_count {
  my ($key, $a, $b) = @_;
  my %d = map {$_ => 1} (split /,/, $a->{'names'}),(split /,/, $b->{'names'});
  my @keys = keys %d;
  return scalar @keys;
}

sub l_join {
  my ($key, $a, $b) = @_;
  return q{} if($a->{$key} eq q{} && $b->{$key} eq q{});
  return $a->{$key} if($a->{$key} ne q{} && $b->{$key} eq q{});
  return $b->{$key} if($a->{$key} eq q{} && $b->{$key} ne q{});
  return $a->{$key}.','.$b->{$key};
}

sub sum_counts {
  my ($key, $a, $b) = @_;
  # ignore $key
  my $sum = 0;
  for(keys %{$b}) {
    next if($_ !~ m/_count$/);
    $sum += $b->{$_};
  }
  return $sum;
}

sub correct_bp {
  my ($key, $a, $b) = @_;
  if($a->{$key} == -1 || $b->{$key} == -1) {
    return -1 if($a->{$key} == $b->{$key});
    die "\n'$key' value combination is impossible for data:\n$a->{_orig}\n$b->{_orig}\n";
  }
  return $b->{'start2'} - $b->{'start1'};
}

sub die_if_diff {
  my ($key, $a, $b) = @_;
  if($a->{$key} ne $b->{$key}) {
    warn "$key, $a->{$key}, $b->{$key}\n";
    warn Dumper($a, $b);
    die "No case provided to handle this difference.\n";
  }
}

sub no_op {
  my ($key, $a, $b) = @_;
  return $b->{$key};
}

sub merge_repeats {
  my ($key, $a, $b) = @_;
  # don't really need the key but use it anyway
  return q{} if($a->{$key} eq q{} && $b->{$key} eq q{});
  return $a->{$key} if($a->{$key} ne q{} && $b->{$key} eq q{});
  return $b->{$key} if($a->{$key} eq q{} && $b->{$key} ne q{});
  return $a->{$key}.','.$b->{$key};
}

sub merge_event {
  my ($this, $that, $header_keys) = @_;
  for my $key(@{$header_keys}) {
    $that->{$key} = $DISPATCH{$key}->($key, $this, $that);
  }
  return 1;
}


1;
