package Sanger::CGP::Brass::Implement;

########## LICENCE ##########
# Copyright (c) 2014 Genome Research Ltd.
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
########## LICENCE ##########


use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use File::Copy qw(copy);
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use Capture::Tiny;
use List::Util qw(first);
use FindBin qw($Bin);
use Capture::Tiny qw(capture_stdout);

use Bio::Brass qw($VERSION);
use File::Copy qw(move);

use PCAP::Threaded;
use PCAP::Bam;

const my $ASSEMBLE_SPLIT => 100;

## input
const my $BAMCOLLATE => q{ outputformat=sam exclude=PROPER_PAIR,UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY mapqthres=6 classes=F,F2 T=%s/bamcollate2_%s filename=%s};
# tmpdir, iter, file
const my $PREP_BAM => q{ -b %s};
# basfile[ -e exclude chrs]
const my $BAMSORT => q{ tmpfile=%s/bamsort_%s inputformat=sam verbose=0 index=1 md5=1 md5filename=%s.md5 indexfilename=%s.bai O=%s};
# out_bamname, out_bamname, out_bamname

## group
const my $BRASS_GROUP => q{ -I %s -F %s %s %s};
# extreme depth, repeats.gz, tumour.brm.bam, normal.brm.bam[, normal_panel.bam]
const my $FILTER_GROUP => q{ -i - -t %s -o %s};
# tumour name
const my $REDIR_GROUP => q{(%s) > %s};
#/software/CGP/projects/brass/bin/brass-group  | /nfs/users/nfs_k/kr2/git/brass/perl/bin/filter_groups.pl -t PD3904a > output.rearr

## filter
const my $BRASS_FILTER => q{ -seq_depth 25.1 -blat %s -ref %s -tumour %s -infile %s -outfile %s};
# path/to/blat, genome.fa, tumour name, groups_in, groups_out, ascat
#perl ~kr2/git/brass/perl/bin/brassI_filter.pl

## assemble
const my $BRASS_ASSEMBLE => q{ -X -m mem -O bedpe -r %s -T %s -o %s %s %s:%s.bai %s:%s.bai};
# genome.fa, tmp, output.tab, groups, tumourbam, tumourbam, normalbam, normalbam

## grass
const my $GRASS => ' -genome_cache %s -ref %s -species %s -assembly %s -platform %s -protocol %s -tumour %s -normal %s -file %s';


sub input {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    ## build command for this index
    #

    my $sample = sanitised_sample_from_bam($input);

    my $outbam = File::Spec->catfile($tmp, sanitised_sample_from_bam($input));
    $outbam .= '.brm.bam';

    my $command = _which('bamcollate2');
    $command .= sprintf $BAMCOLLATE, $tmp, $index, $input;
    $command .= ' | ';
    $command .= "$^X ";
    $command .= _which('brassI_prep_bam.pl');
    $command .= sprintf $PREP_BAM, $input.'.bas';
    $command .= " -e $options->{'exclude'}" if(exists $options->{'exclude'});
    $command .= ' | ';
    $command .= _which('bamsort');
    $command .= sprintf $BAMSORT, $tmp, $index, $outbam, $outbam, $outbam;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub group {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups');
  my $tumour = File::Spec->catfile($tmp, $options->{'safe_tumour_name'}).'.brm.bam';
  my $normal = File::Spec->catfile($tmp, $options->{'safe_normal_name'}).'.brm.bam';

  my $command = _which('brass-group');
  $command .= sprintf $BRASS_GROUP, $options->{'depth'}, $options->{'repeats'}, $tumour, $normal;
  $command .= ' | ';
  $command .= "$^X ";
  $command .= _which('brassI_pre_filter.pl');
  $command .= sprintf $FILTER_GROUP, $options->{'tumour_name'}, $groups;
  $command .= " -n $options->{filter}" if(exists $options->{'filter'});

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub filter {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $decomp = File::Spec->catfile($tmp, 'genome.fa');
  unless(-e $decomp) {
    if($options->{genome} =~ m/\.[gz|razf]$/) {
      system([0,2], "gunzip -c $options->{genome} > $decomp");
      system("samtools faidx $decomp");
    }
    else {
      symlink($options->{genome}, $decomp);
      symlink("$options->{genome}.fai", "$decomp.fai");
    }
  }

  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups');
  my $filtered = $groups.'.filtered';
  my $blat = _which('blat');
  my $command = "$^X ";
  $command .= _which('brassI_filter.pl');
  $command .= sprintf $BRASS_FILTER, $blat,
                                      $decomp,
                                      $options->{'safe_tumour_name'},
                                      $groups,
                                      $filtered;
  $command .= " -ascat $options->{ascat}" if(exists $options->{'ascat'} && defined $options->{'ascat'});

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub split_filtered {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);


  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups.filtered.bedpe');

  my $split_dir = File::Spec->catdir($tmp, 'split');
  remove_tree($split_dir) if(-d $split_dir);
  make_path($split_dir);
  my $command = sprintf 'split --suffix-length=3 --numeric-suffixes --verbose --lines=%s %s %s/split.',
                        $ASSEMBLE_SPLIT, $groups, $split_dir;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub assemble {
 my ($index_in, $options) = @_;

  my $tmp = $options->{'tmp'};

	# first handle the easy bit, skip if limit not set
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});


  my @indicies = limited_indices($options, $index_in, $options->{'splits'});
  for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

    my $split_dir = File::Spec->catdir($tmp, 'split');
    my $split_file = File::Spec->catfile($split_dir, 'split.');
    $split_file .= sprintf '%03d', $index-1;

    my $tmp_assemble = File::Spec->catdir($tmp, 'assemble');
    make_path($tmp_assemble) unless(-e $tmp_assemble);

    my $assembled = File::Spec->catfile($tmp_assemble, 'bedpe.');
    $assembled .= sprintf '%03d', $index-1;

    my $command = "$^X ";
    $command .= _which('brass-assemble');
    $command .= sprintf $BRASS_ASSEMBLE, $options->{'genome'},
                                          $tmp_assemble,
                                          $assembled,
                                          $split_file,
                                          $options->{'tumour'}, $options->{'tumour'},
                                          $options->{'normal'}, $options->{'normal'};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub grass {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $assembled = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.assembled.bedpe');
  my $merge = sprintf '(cat %s/assemble/bedpe.* | sort -k1,1 -k 2,2n > %s)', $tmp, $assembled;

  my $tumour = File::Spec->catfile($tmp, $options->{'safe_tumour_name'}).'.brm.bam';
  my $normal = File::Spec->catfile($tmp, $options->{'safe_normal_name'}).'.brm.bam';

  my $grass_cmd = "$^X ";
  $grass_cmd .= _which('grass.pl');
  $grass_cmd .= sprintf $GRASS, $options->{'g_cache'},
                              $options->{'genome'},
                              $options->{'species'},
                              $options->{'assembly'},
                              $options->{'platform'},
                              $options->{'protocol'},
                              $options->{'tumour_name'},
                              $options->{'normal_name'},
                              $assembled;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), [$merge, $grass_cmd], 0);

  my $munged_name =File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'_ann.assembled.vcf');
  my $annotated = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.annot.vcf');
  move $munged_name, $annotated || die $!;
  $munged_name =~ s/vcf$/bedpe/;
  $annotated =~ s/vcf$/bedpe/;
  move $munged_name, $annotated || die $!;

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return 1;
}

sub split_count {
  my $options = shift;
  my $split_dir = File::Spec->catdir($options->{'tmp'}, 'split');
  my $split_count = 0;
  opendir(my $dh, $split_dir);
  while(readdir $dh) {
    $split_count++ if($_ =~ m/^split\.[[:digit:]]+$/);
  }
  closedir $dh;
  return $split_count;
}

sub limited_split {
  my ($options, $index_in) = @_;
	my $split_count = split_count($options);
	return limited_indices($options, $index_in, $split_count);
}

sub limited_indices {
	my ($options, $index_in, $count) = @_;
  my @indicies;
  if(exists $options->{'limit'}) {
    # main script checks index is not greater than limit or < 1
	  my $base = $index_in;
	  while($base <= $count) {
	    push @indicies, $base;
	    $base += $options->{'limit'};
	  }
	}
	else {
	  push @indicies, $index_in;
	}
	return @indicies;
}

sub tabix {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $annotated = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.annot.vcf');
  my $sorted = $annotated.'.srt';

  my $header = sprintf q{(grep '^#' %s > %s)}, $annotated, $sorted;
  my $sort = sprintf q{(grep -v '^#' %s | sort -k 1,1 -k 2,2n >> %s)}, $annotated, $sorted;

  my $vcf_gz = $annotated.'.gz';
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $sorted, $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  my @commands = ($header, $sort, $bgzip, $tabix);

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return 1;
}

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  $options->{'safe_tumour_name'} = sanitised_sample_from_bam($options->{'tumour'});
  $options->{'safe_normal_name'} = sanitised_sample_from_bam($options->{'normal'});
  return 1;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
  return $path;
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_-]/_/ig; # sanitise sample name
  return $sample;
}

sub get_bam_info {
  my $options = shift;
  my $bam_ob = PCAP::Bam->new($options->{'tumour'});

  my ($tum_name, $tum_plat);
  for(@{$bam_ob->read_group_info()}) {
    if(exists $_->{'SM'}) {
      $tum_name = $_->{'SM'} unless(defined $tum_name);
      die "Tumour BAM file has multiple sample names in header ($options->{tumour})\n" if($tum_name ne $_->{'SM'});
    }
    if(exists $_->{'PL'}) {
      $tum_plat = $_->{'PL'} unless(defined $tum_plat);
      die "Tumour BAM file has multiple sequencing platforms in header ($options->{tumour})\n" if($tum_plat ne $_->{'PL'});
    }
  }

  $bam_ob = PCAP::Bam->new($options->{'normal'});
  my ($norm_name, $norm_plat);
  for(@{$bam_ob->read_group_info()}) {
    if(exists $_->{'SM'}) {
      $norm_name = $_->{'SM'} unless(defined $norm_name);
      die "Normal BAM file has multiple sample names in header ($options->{normal})\n" if($norm_name ne $_->{'SM'});
    }
    if(exists $_->{'PL'}) {
      $norm_plat = $_->{'PL'} unless(defined $norm_plat);
      die "Normal BAM file has multiple sequencing platforms in header ($options->{normal})\n" if($norm_plat ne $_->{'PL'});
    }
  }

  my $platform;
  if(defined $tum_plat && defined $norm_plat) {
    die "BAMs have different sequencing platform: $tum_plat, $norm_plat\n" if($norm_plat ne $tum_plat);
    $platform = $norm_plat;
  }

  die "Specified tumour name doesn't match BAM header" if(defined $options->{'tumour_name'} && defined $tum_name && $options->{'tumour_name'} ne $tum_name);
  die "Specified normal name doesn't match BAM header" if(defined $options->{'normal_name'} && defined $norm_name && $options->{'normal_name'} ne $norm_name);

  $options->{'tumour_name'} = $tum_name unless(defined $options->{'tumour_name'});
  $options->{'normal_name'} = $norm_name unless(defined $options->{'normal_name'});

  die "Specified platform name doesn't match BAM headers" if(defined $options->{'platform'} && defined $platform && $options->{'platform'} ne $platform);
  $options->{'platform'} = $platform unless(defined $options->{'platform'});

  return $options;
}

1;
