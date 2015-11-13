#!/usr/bin/perl

use strict;

while(<>) {
  my ($flag, $this_chr, $mate_chr) = (split /\t/)[1,2,6];
  if(($mate_chr eq q{=} || $mate_chr eq $this_chr)
    &&
    ((($flag & 16) && ($flag & 32)) || (!($flag & 16) && !($flag & 32)))
    ) { print; }
}
