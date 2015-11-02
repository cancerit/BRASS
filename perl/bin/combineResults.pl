#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use autodie;
use Const::Fast qw(const);

const my @NEW_BEDPE_HEADER => ('# chr1','start1','end1','chr2','start2','end2','id/name','brass_score','assembly_score','strand1','strand2','sample','readpair names','readpair count','bal_trans','inv','occL','occH','copynumber_flag','range_blat','Brass Notation','non-template','micro-homology','assembled readnames','assembled read count','gene','gene_id','transcript_id','strand','end_phase','region','region_number','total_region_count','first/last','gene','gene_id','transcript_id','strand','phase','region','region_number','total_region_count','first/last','fusion_flag');

if(scalar @ARGV != 3) {
  warn "USAGE: combineResults.pl X_ann.groups X_ann.assembled X.final\n";
  warn "\tX.final is a prefix, relevant suffixes will be added for VCF and BEDPE outputs\n";
  exit 1;
}
my ($groups_prefix, $assembled_prefix, $final_prefix) = @ARGV;


mergeVcf($groups_prefix, $assembled_prefix, $final_prefix);
mergeBedpe($groups_prefix, $assembled_prefix, $final_prefix);

sub mergeVcf {
  my ($groups_prefix, $assembled_prefix, $final_prefix) = @_;

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

  open my $FINAL, '>', $final_vcf;
  open my $FIXED, '<', $phaseI_vcf;

  my @ends;

  while(my $line = <$FIXED>) {
    if(@ends == 2) {
      svclass(\@ends);
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
      my $old_info = $orig_data{$id}->[7];
      my ($tsrds) = $bits[7] =~ m/(TSRDS=[^;]+)/;

      $bits[7] = $old_info.';'.$tsrds;
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
    svclass(\@ends);
    print $FINAL join("\t",@{$ends[0]}), "\n";
    print $FINAL join("\t",@{$ends[1]}), "\n";
  }

  close $FIXED;
  close $FINAL;
}

sub mergeBedpe {
  my ($groups_prefix, $assembled_prefix, $final_prefix) = @_;
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
      die "merged events correlate with assembled results, don't know what to do";
    }
    elsif(!exists $orig_data{$id}) {
      my @new;
      push @new, @bits[0..7]; # 1-8
      push @new, q{_}; # 9
      push @new, @bits[8..9]; # 10-11
      push @new, $bits[16]; # 12
      push @new, @bits[18..25]; # 13-20
      push @new, q{_},q{_},q{_},q{_},q{_}; # 21-25
      push @new, @bits[26..44]; # 26-44

      print $FINAL join("\t", @new),"\n";
      next;
    }
    else {
      my @old_brass_II = @{$orig_data{$id}};
      my @new;
      push @new, @old_brass_II[0..6]; # 1-7
      push @new, q{_}; # 8
      push @new, @old_brass_II[7..10]; # 9-12
      push @new, @bits[18..25]; # 13-20
      push @new, @old_brass_II[11..14]; # 21-24
      push @new, scalar (split /,/, $new[-1]); # 25 # assembled read count
      push @new, @bits[26..44]; # 26-44

      print $FINAL join("\t", @new),"\n";
      next;
    }
  }


  close $FIXED;
  close $FINAL;

  return;
}


sub svclass {
  my ($end_a, $end_b) = @{$_[0]};
  die "Record IDs out of sync at:\n",join("\t",@{$end_a}),"\n",join("\t",@{$end_b}),"\n" if($end_a->[2] !~ m/_[12]$/ || $end_b->[2] !~ m/_[12]$/);
  my $class;
  if($end_a->[0] ne $end_b->[0]) {
    $class = 'translocation';
  }
  elsif($end_a->[4] =~ m/^[[:upper:]]\[/) {
    # ++
    $class = 'deletion';
  }
  elsif($end_a->[4] =~ m/^[[:upper:]]\]/ || $end_a->[4] =~ m/\[[[:upper:]]$/) {
    # +- / -+
    $class = 'inversion';
  }
  elsif($end_a->[4] =~ m/\][[:upper:]]$/) {
    # --
    $class = 'tandem-duplication';
  }
  else {
    die "Unknown rearrangement syntax: $end_a->[4]\n";
  }
  $end_a->[7] .= ';SVCLASS='.$class;
  $end_b->[7] .= ';SVCLASS='.$class;
  return;
}

