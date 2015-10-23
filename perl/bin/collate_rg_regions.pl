#!/usr/bin/perl

use strict;
use List::Util qw(first);

my ($full_bed_pe, $rg_patt_pe, $out_file) = @ARGV;

my @ids;
open my $FILT, '<', $rg_patt_pe || die $!;
while(my $line = <$FILT>) {
  push @ids, (split /\t/, $line)[6];
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
  my $id = (split /\t/, $line)[6];
  print $line if(first{$id == $_} @ids);
}
close $MAIN;

close $ofh if(defined $out_file);
