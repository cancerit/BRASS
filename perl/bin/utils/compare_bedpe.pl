#!/usr/bin/perl

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

use strict;
use Data::Dumper;
use IO::Uncompress::Gunzip qw($GunzipError);
use List::Util qw(first);
use Data::Dumper;

if(@ARGV < 2) {
  warn "USAGE: ./compare_bedpe.pl old.annot.bedpe.gz new.annot.bedpe.gz [new.r6 <CLASS_STR>]\n\n";
  warn "\t'*.r6' will be in the intermediates.tar.gz\n";
  warn "\tCLASS_STR is constructed by concatenating 'Type' and 'Support' with ':' to list out the IDs from that grouping of results.\n";
  warn "\t\t- run without this set to see possible values.\n\n";
  exit 1;
}

# old, new
my ($file_old, $file_new, $r6, $dump_class) = @ARGV;

my $plus_split_read_info = 0;
if(defined $r6) {
  $plus_split_read_info = 1;
}


my $grps_old = load_groups($file_old, 1);
my $grps_new = load_groups($file_new);
my $new_coord = load_coord($r6);

compare_grps($grps_old, $grps_new, $new_coord);

sub compare_grps {
  my ($old, $new, $coord) = @_;

  my %all_criteria;

  my ($common, $common_pass, $common_split, $common_ass, $new_only, $new_pass, $new_split, $new_ass, $multi_hit) = (0,0,0,0,0,0);

  # only compare events from new, if can't find overlap there may be a problem
  for my $chr_str(sort keys %{$new}) {
    for my $region_n(sort {$a->[0] <=> $b->[0]} @{$new->{$chr_str}}) {
      my $cn_pass = $region_n->[8];
      my ($l_start_new, $l_end_new, $h_start_new, $h_end_new, $brass_new, $assembly_new, $id_new, $svc_new) = @{$region_n};
      push @{$all_criteria{$id_new}}, 'CN' if($cn_pass == 1);
      push @{$all_criteria{$id_new}}, 'ASSEMBLE' if($assembly_new ne q{_});
      printf "\n%-7s: %d-%d : %d - %d %s (%2s/%-3s)", $chr_str, $l_start_new, $l_end_new, $h_start_new, $h_end_new, $svc_new, $brass_new, $assembly_new;
      if($plus_split_read_info && exists $coord->{$id_new}) {
        printf "\t(split_read info: %s)\n", join(q{,}, @{$coord->{$id_new}});
        push @{$all_criteria{$id_new}}, 'SPLIT';
      }
      else {
        print "\n";
      }
      my $old_found = 0;
      if(exists $old->{$chr_str}) {
        for my $region_o(sort {$a->[0] <=> $b->[0]}  @{$old->{$chr_str}}) {
          my ($l_start_old, $l_end_old, $h_start_old, $h_end_old, $brass_old, $assembly_old, $id_old, $svc_old) = @{$region_o};
          next if($l_start_old > $l_end_new+500);
          next if($l_end_old < $l_start_new-500);
          next if($h_start_old > $h_end_new+500);
          next if($h_end_old < $h_start_new-500);
          printf "%-7s: %d-%d : %d - %d %s (%2s/%-3s)\n", q{Old}, $l_start_old, $l_end_old, $h_start_old, $h_end_old, $svc_old, $brass_old, $assembly_old;
          $old_found++;
        }
      }
      if($old_found == 0) {
        push @{$all_criteria{$id_new}}, 'UNIQUE';
        $new_only++;
        $new_pass++ if($cn_pass == 1);
      }
      elsif($old_found == 1) {
        push @{$all_criteria{$id_new}}, 'COMMON';
        $common++;
        $common_pass++ if($cn_pass == 1);
      }
      else {
        $multi_hit++;
      }
    }
  }

  my $grouped = classify(\%all_criteria);
  if(defined $dump_class) {
    print 'IDs: ';
    print join q{,}, sort {$a<=>$b} @{$grouped->{$dump_class}};
    print "\n";
  }
}

sub classify {
  my $result = shift;
  my @ids = keys %{$result};
  my %stats;

  my $common = 0;
  my $unique = 0;
  my $total = 0;

  my %final_groups;

  for my $id(@ids) {
    my $tag = join(q{,}, sort @{$result->{$id}}, 'GROUPS');

    if($tag =~ s/COMMON//) {
      $tag =~ s/^,//;
      $tag =~ s/,$//;
      $tag =~ s/,{2}/,/g;
      $stats{'COMMON'}{$tag}++;
      $common++;
      $total++;

      push @{$final_groups{'COMMON:'.$tag}}, $id;
    }
    else {
      $tag =~ s/UNIQUE//;
      $tag =~ s/^,//;
      $tag =~ s/,$//;
      $tag =~ s/,{2}/,/g;
      $stats{'UNIQUE'}{$tag}++;
      $unique++;
      $total++;

      push @{$final_groups{'UNIQUE:'.$tag}}, $id;
    }



  }

  print "Overall:\n";
  print "\tTotal: $total\n";
  print "\tCommon: $common\n";
  print "\tUnique: $unique\n";

  print "Detail:\n\t#Count\tType\tSupport\n";
  for my $top(sort keys %stats) {
    for my $mid (sort keys %{$stats{$top}}) {
      print "\t$stats{$top}{$mid}\t$top\t$mid\n";
    }
  }
  return \%final_groups;
}

sub load_groups {
  my ($file, $fix_ccords) = @_;
  my %groups;
  my $z = IO::Uncompress::Gunzip->new($file, MultiStream => 1)
            or die "gunzip failed: $GunzipError\n";
  while(my $line = <$z>) {
    next if($line =~ m/^#/);
    chomp $line;
    my @tmp = split /\t/, $line;
    my $cn_flag = $tmp[20];
    my ($l_c, $l_s, $l_e, $h_c, $h_s, $h_e, $id, $score_brass, $l_d, $h_d, $samples, $svclass, undef, $score_ass, @the_rest) = @tmp;

next if($samples =~ m/,/);

    if($l_d eq q{-}) {
      $l_s++;
      ($l_s, $l_e) = ($l_e, $l_s);
      $l_s--;
    }
    if($h_d eq q{-}) {
      $h_s++;
      ($h_s, $h_e) = ($h_e, $h_s);
      $h_s--;
    }
    my $toplev = sprintf '%s%s/%s%s', $l_c, $l_d, $h_c, $h_d;
    push @{$groups{$toplev}}, [$l_s, $l_e, $h_s, $h_e, $score_brass, $score_ass, $id, $svclass, $cn_flag, \@the_rest];
  }
  close $z;
  return \%groups;
}

sub load_coord {
  my ($file) = @_;
  my %groups;
  open my $fh, '<', $file || die $!;
  while(my $line = <$fh>) {
    next if($line =~ m/^#/);
    chomp $line;
    my ($id, $low, $high, @the_rest) = (split /\t/, $line)[6,12,13];#,14,15,16,17];
    $groups{$id} = [$low, $high, @the_rest];
  }
  close $fh;
  return \%groups;
}
