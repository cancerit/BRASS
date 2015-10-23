#!/usr/bin/perl

use strict;
use autodie;
use warnings FATAL => 'all';
use Cwd;
use Try::Tiny qw(try finally);
use Const::Fast qw(const);

const my $FILE_PAIR => '%s.%s.ngscn.bed.gz';
const my $FILE_FOLD => '%s.%s.ngscn.fb_reads.bed.gz';

die "Usage: genome.fa.fai sample_name indir" unless(scalar @ARGV == 3);

my ($fai, $sample, $indir) = @ARGV;
die "ERROR: *.fai file must exist with non-zero size\n" unless(-e $fai && -s _);
die "ERROR: indir must exist\n" unless(-e $indir && -d _);

my @chr_order;
open my $FAI_IN, '<', $fai;
while(<$FAI_IN>) {
  push @chr_order, (split /\t/)[0];
}
close $FAI_IN;

my $final_pair = "$indir/$sample.ngscn.bed.gz";
my $final_fold = "$indir/$sample.ngscn.fb_reads.bed.gz";

unlink $final_pair if(-e $final_pair);
unlink $final_fold if(-e $final_fold);

my $init_dir = getcwd;

try {
  chdir $indir;
  zcat_to_gzip($FILE_PAIR, $final_pair, $sample, \@chr_order);
  zcat_to_gzip($FILE_FOLD, $final_fold, $sample, \@chr_order);
} finally {
  chdir $init_dir;
};

sub zcat_to_gzip {
  my ($format, $outfile, $sample, $chrs) = @_;
  my @args;
  for my $chr(@{$chrs}) {
    push @args, sprintf $format, $sample, $chr;
    die "Expected file missing $indir/$args[-1]\n" unless(-e $args[-1]);
  }
  system("bash -c 'set -o pipefail; zcat @args | gzip -c > $outfile'") == 0 or die "Failed to merge files to $outfile: $!\n";
}
