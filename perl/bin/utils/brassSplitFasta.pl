#!/usr/bin/perl

use strict;

use constant FILES => 4;

my ($fasta, $outdir) = @ARGV;

$outdir =~ s|/$||;

warn "Calculating records per file...\n";

my $records = `grep -c '^>' $fasta`;
chomp $records;

my $records_per_file = int (($records / FILES) + 1);

warn "Total records: $records\n";
warn "Records per file: $records_per_file (".FILES.")\n";

$records = 0;
my $file_count = 1;

my $OUT;

open my $IN, '<', $fasta or $!;
while(my $line = <$IN>) {
  if($line =~ m/^>/) {
    $records = 0 if($records >= $records_per_file);
    if($records == 0) {
      close $OUT if(defined $OUT);
      my $outfile = "$outdir/all_ncbi_bacteria.${file_count}.fa";
      warn "Generating: $outfile\n";
      open $OUT, '>', $outfile or $!;
      $file_count++;
    }
    $records++;
  }
  print $OUT $line;
}
close $IN;
close $OUT;
