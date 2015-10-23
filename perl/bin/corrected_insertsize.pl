#!/usr/bin/perl

use strict;

my @count_of;
while(<>) {
  my ($this_chr, $mate_chr, $isize) = (split /\t/)[2,6,8];
  if(($mate_chr eq q{=} || $mate_chr eq $this_chr)
    &&
    $isize >= 0
    ) { $count_of[$isize]++; }
}

for(0..3000) {
  print "$_ $count_of[$_]\n";
}
