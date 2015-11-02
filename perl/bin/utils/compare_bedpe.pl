#!/usr/bin/perl

use strict;
use Data::Dumper;

# old, new
my ($file_old, $file_new, $r6) = @ARGV;

my $grps_old = load_groups($file_old, 1);
my $grps_new = load_groups($file_new);
my $new_coord = load_coord($r6);

compare_grps($grps_old, $grps_new, $new_coord);

sub compare_grps {
  my ($old, $new, $coord) = @_;
  # only compare events from new, if can't find overlap there may be a problem
  for my $chr_str(sort keys %{$new}) {
    next unless(exists $old->{$chr_str});
    for my $region_n(sort {$a->[0] <=> $b->[0]} @{$new->{$chr_str}}) {
      my ($l_start_new, $l_end_new, $h_start_new, $h_end_new, $brass_new, $assembly_new, $id_new) = @{$region_n};
      print sprintf "\n%-7s: %d-%d : %d - %d (%2s/%-3s)", $chr_str, $l_start_new, $l_end_new, $h_start_new, $h_end_new, $brass_new, $assembly_new;
      print sprintf "\t(split_read info: %s)\n", join(q{,}, @{$coord->{$id_new}});
#      print sprintf "\t!! NEW low coords reversed !!\n" if($l_start_new > $l_end_new);
#      print sprintf "\t!! NEW high coords reversed !!\n" if($h_start_new > $h_end_new);
      for my $region_o(sort {$a->[0] <=> $b->[0]}  @{$old->{$chr_str}}) {
        my ($l_start_old, $l_end_old, $h_start_old, $h_end_old, $brass_old, $assembly_old, $id_old) = @{$region_o};
        next if($l_start_old > $l_end_new+500);
        next if($l_end_old < $l_start_new-500);
        print sprintf "%-8s %d-%d : %d - %d (%2s/%-3s)\n", q{Old}, $l_start_old, $l_end_old, $h_start_old, $h_end_old, $brass_old, $assembly_old;
#        print sprintf "\t\t!! OLD low coords reversed !!\n" if($l_start_old > $l_end_old);
#        print sprintf "\t\t!! OLD high coords reversed !!\n" if($h_start_old > $h_end_old);
      }
    }
  }
}

sub load_groups {
  my ($file, $fix_ccords) = @_;
  my %groups;
  open my $fh, '<', $file || die $!;
  while(my $line = <$fh>) {
    next if($line =~ m/^#/);
    chomp $line;
    my ($l_c, $l_s, $l_e, $h_c, $h_s, $h_e, $id, $score_brass, $score_ass, $l_d, $h_d, @the_rest) = split /\t/, $line;
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
    push @{$groups{$toplev}}, [$l_s, $l_e, $h_s, $h_e, $score_brass, $score_ass, $id, \@the_rest];
  }
  close $fh;
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
