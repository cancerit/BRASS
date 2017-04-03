#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
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
use warnings FATAL => 'all';
use autodie;
use List::Util qw(first);
use Const::Fast qw(const);

const my @CN_STATES => qw(FAIL-DEFAULT DEFAULT REAL);
const my @NEW_BEDPE_HEADER => ('# chr1','start1','end1','chr2','start2','end2','id/name','brass_score','strand1','strand2','sample','svclass','bkdist','assembly_score','readpair names','readpair count','bal_trans','inv','occL','occH','copynumber_flag','range_blat','Brass Notation','non-template','micro-homology','assembled readnames','assembled read count','gene1','gene_id1','transcript_id1','strand1','end_phase1','region1','region_number1','total_region_count1','first/last1','gene2','gene_id2','transcript_id2','strand2','phase2','region2','region_number2','total_region_count2','first/last2','fusion_flag');

if(scalar @ARGV < 3) {
  warn "USAGE: combineResults.pl X_ann.groups X_ann.assembled X.final ploidy Acf cn_state\n";
  warn "\tX.final is a prefix, relevant suffixes will be added for VCF and BEDPE outputs\n";
  exit 1;
}
my ($groups_prefix, $assembled_prefix, $final_prefix, $ploidy, $acf, $cn_state) = @ARGV;

die "ERROR: cn_state must be one of ".join(',',@CN_STATES)." - you provided $cn_state\n" unless(first { $cn_state eq $_ } @CN_STATES);

my $ascat_info = sprintf "Ploidy=%s,Acf=%s,CnState=%s", $ploidy, $acf, $cn_state;
my $svclass_bkpt_dist = mergeBedpe($groups_prefix, $assembled_prefix, $final_prefix, $ascat_info);
mergeVcf($groups_prefix, $assembled_prefix, $final_prefix, $ascat_info, $svclass_bkpt_dist);


sub mergeVcf {
  my ($groups_prefix, $assembled_prefix, $final_prefix, $ascat_info, $svclass_bkpt_dist) = @_;

  my $assembled_vcf = "$assembled_prefix.vcf";
  my $phaseI_vcf = "$groups_prefix.vcf";

  my $final_vcf = "$final_prefix.vcf";

  my %orig_data;
  open my $ORIG, '<', $assembled_vcf;
  while(my $line = <$ORIG>) {
    next if($line =~ m/^#/);
    chomp $line;
    my @bits = split /\t/, $line;
    my $id = $bits[2];
    $orig_data{$id} = \@bits;
  }
  close $ORIG;

  my $new_info = q{##INFO=<ID=BAS,Number=.,Type=Integer,Description="Brass Assembly Score:A maximum score of 100 indicates a perfect pattern of 5 vertices in Velvet's de Bruijn graph">}
                  .qq{\n}
                  .q{##INFO=<ID=SVCLASS,Number=.,Type=String,Description="basic SV class, deletion, inversion, tandem-duplication, translocation">};
  my $new_format = q{##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Count of pairs that span this breakpoint">};

  my $info_found = 0;
  my $format_found = 0;
  my $sample_found = 0;

  open my $FINAL, '>', $final_vcf;
  open my $FIXED, '<', $phaseI_vcf;

  my @ends;

  while(my $line = <$FIXED>) {
    if(@ends == 2) {
      my $id = $ends[0]->[2];
      my ($id_base) = (split q{_}, $id)[0];
      for my $end(@ends) {
        $end->[7] =~ s/BKDIST=[[:digit:]]+//;
        $end->[7] .= ';BKDIST='.$svclass_bkpt_dist->{$id_base}->{'bkptdist'};
        $end->[7] .= ';SVCLASS='.$svclass_bkpt_dist->{$id_base}->{'svclass'};
        $end->[7] =~ s/;{2}/;/g;
      }
      print $FINAL join("\t",@{$ends[0]}), "\n";
      print $FINAL join("\t",@{$ends[1]}), "\n";
      @ends = ();
    }
    if($line =~ m/^#/) {
      if($info_found == 0 && $line =~ m/^##INFO/) {
        print $FINAL $new_info,"\n";
        $info_found = 1;
      }
      if($format_found == 0 && $line =~ m/^##FORMAT/) {
        print $FINAL $new_format,"\n";
        $format_found = 1;
      }
      if($sample_found == 0 && $line =~ m/^##SAMPLE.*=TUMOUR,/) {
        chomp $line;
        my $trailer = chop $line;
        $line = sprintf "%s,%s%s\n", $line, $ascat_info, $trailer;
        $sample_found = 1;
      }
      print $FINAL $line;
      next;
    }
    chomp $line;
    my @bits = split /\t/, $line;
    my $id = $bits[2];

    if(!exists $orig_data{$id}) {
      $bits[-3] .= ':PS';
      $bits[-2] = '0:'.$bits[-2]; # normal is always 0:0
      $bits[-1] = '0:'.$bits[-1];

      push @ends, \@bits;
      next;
    }
    else {
      my $original_info = $orig_data{$id}->[7];
      my ($tsrds) = $bits[7] =~ m/(TSRDS=[^;]+)/;
      $bits[1] = $orig_data{$id}->[1];
      $bits[3] = $orig_data{$id}->[3];
      $bits[4] = $orig_data{$id}->[4];
      $bits[7] = $original_info.';'.$tsrds;
      $bits[7] .= ';BAS='.$orig_data{$id}->[5];
      $bits[5] = q{.};
      $bits[-3] .= ':PS';
      $bits[-2] = '0:'.$bits[-2]; # normal is always 0:0
      $bits[-1] = $orig_data{$id}->[-1].':'.$bits[-1];

      $bits[7] =~ s/;;/;/g;

      push @ends, \@bits;
      next;
    }
  }
  if(@ends == 2) {
    my $id = $ends[0]->[2];
    my ($id_base) = (split q{_}, $id)[0];
    for my $end(@ends) {
      $end->[7] =~ s/BKDIST=[[:digit:]]+//;
      $end->[7] .= ';BKDIST='.$svclass_bkpt_dist->{$id_base}->{'bkptdist'};
      $end->[7] .= ';SVCLASS='.$svclass_bkpt_dist->{$id_base}->{'svclass'};
      $end->[7] =~ s/;{2}/;/g;
    }
    print $FINAL join("\t",@{$ends[0]}), "\n";
    print $FINAL join("\t",@{$ends[1]}), "\n";
  }

  close $FIXED;
  close $FINAL;
}

sub mergeBedpe {
  my ($groups_prefix, $assembled_prefix, $final_prefix, $ascat_info) = @_;
  my %svclass_bkpt_dist;
  my $assembled_bedpe = "$assembled_prefix.bedpe";
  my $phaseI_bedpe= "$groups_prefix.bedpe";
  my $final_bedpe = "$final_prefix.bedpe";

  my %orig_data;
  open my $ORIG, '<', $assembled_bedpe;
  while(my $line = <$ORIG>) {
    next if($line =~ m/^#/);
    chomp $line;
    my @bits = split /\t/, $line;
    my $id = $bits[6];
    $orig_data{$id} = \@bits;
  }
  close $ORIG;

  open my $FINAL, '>', $final_bedpe;
  open my $FIXED, '<', $phaseI_bedpe;

  print $FINAL sprintf "##%s\n", $ascat_info;

  while(my $line = <$FIXED>) {
    if($line =~ m/^#/) {
      if($line =~ m/^# chr/) {
        print $FINAL join("\t", @NEW_BEDPE_HEADER),"\n";
      }
      else {
        print $FINAL $line;
      }
      next;
    }
    chomp $line;
    my @bits = split /\t/, $line;
    my $id = $bits[6];
    my @new;
    if(index($id, ',') != -1) {
      my @ids = split /,/, $id;
      my @multi_record;
      for(@ids) {
        push @multi_record, $orig_data{$id} if(exists $orig_data{$id});
      }
      if(scalar @multi_record) {
        warn "MULTI\n";
        warn "$line\n";
        warn Dumper(\@multi_record);
      }
      die q{merged events correlate with assembled results, don't know what to do};
    }
    elsif(!exists $orig_data{$id}) {
      # did not assemble
      push @new, @bits[0..7]; # 1-8
      push @new, @bits[8..9]; # 9-10
      push @new, $bits[16]; # 11
      push @new, svclass_bedpd(@new[0,3,8,9]); # 12
      push @new, svdist_bedpd(@new[0..5]); # 13
      push @new, q{_}; # 14 assembly_score
      push @new, @bits[18..25]; # 15-22
      push @new, q{_},q{_},q{_},q{_},q{_}; # 23-27
      push @new, @bits[26..44]; # 28-46
    }
    else {
      my @old_brass_II = @{$orig_data{$id}};
      push @new, @old_brass_II[0..6]; # 1-7
      # original brass_score here!!
      push @new, $bits[7]; # 8
      push @new, @old_brass_II[8..10]; # 9-11
      push @new, svclass_bedpd(@new[0,3,8,9]); # 12
      push @new, svdist_bedpd(@new[0..5]); # 13
      push @new, $old_brass_II[7]; # 14 assembly_score
      push @new, @bits[18..25]; # 15-22
      push @new, @old_brass_II[11..14]; # 23-26
      push @new, scalar (split /,/, $new[-1]); # 27 # assembled read count
      push @new, @bits[26..44]; # 28-46
    }
    clean_record(\@new);
    print $FINAL join("\t", @new),"\n";
    $svclass_bkpt_dist{$id} = { 'svclass' => $new[11],
                                'bkptdist' =>$new[12],};
  }

  close $FIXED;
  close $FINAL;

  return \%svclass_bkpt_dist;
}

sub clean_record {
  my $rec = shift;
  for(0..@{$rec}-1) {
      $rec->[$_] = '_' unless(defined $rec->[$_]);
  }
}

sub svclass_bedpd {
  my ($chrL, $chrH, $strL, $strH) = @_;
  return 'translocation' if($chrL ne $chrH);
  return 'inversion' if($strL ne $strH);
  return 'deletion' if($strL eq '+' && $strH eq '+');
  return 'tandem-duplication' if($strL eq '-' && $strH eq '-');
}

sub svdist_bedpd {
  my ($chrL, $startL, $endL, $chrH, $startH, $endH) = @_;
  return -1 if($chrL ne $chrH);
  return abs $startH - $endL;
}
