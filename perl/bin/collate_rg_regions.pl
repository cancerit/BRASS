#!/usr/bin/perl

use strict;

my ($full_bed_pe, $rg_patt_pe, $out_file) = @ARGV;

my %ids;
open my $FILT, '<', $rg_patt_pe || die $!;
while(my $line = <$FILT>) {
  my ($id, $low, $high) = (split /\t/, $line)[6,12,13];
  $ids{$id} = [$low, $high];
}
close $FILT;

my $ofh = *STDOUT;
if(defined $out_file) {
  open $ofh, '>', $out_file || die $!;
}

open my $MAIN, '<', $full_bed_pe || die $!;
while(my $line = <$MAIN>) {
  if($line =~ m/^#/) {
    print $line;
    next;
  }

  my @F = split /\t/, $line;
  my $id = $F[6];
  next unless(exists $ids{$id});
  $F[1] = $ids{$id}->[0]-1;
  $F[2] = $ids{$id}->[0];
  $F[4] = $ids{$id}->[1]-1;
  $F[5] = $ids{$id}->[1];
  print join("\t", @F);
}
close $MAIN;

close $ofh if(defined $out_file);
