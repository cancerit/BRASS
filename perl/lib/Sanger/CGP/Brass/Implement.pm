package Sanger::CGP::Brass::Implement;

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use autodie qw(:all);
use Const::Fast qw(const);
use File::Copy qw(copy move);
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile tempdir);
use Capture::Tiny;
use List::Util qw(first);
use FindBin qw($Bin);
use Capture::Tiny qw(capture_stdout capture);
use File::ShareDir qw(module_dir);
use Cwd qw(abs_path getcwd);
use File::Basename;

use Bio::Brass qw($VERSION);
use File::Copy qw(move);

use PCAP::Threaded;
use PCAP::Bam;

const my $ASSEMBLE_SPLIT => 40;

## input
const my $BED_INTERSECT_V => q{ intersect -ubam -v -abam %s -b %s };
const my $BAMCOLLATE => q{ outputformat=sam exclude= classes=F,F2 T=%s/bamcollate2_%s};
# tmpdir, iter, file
const my $PREP_BAM => q{ -b %s -p %s};
# basfile[ -e exclude chrs]
const my $BAMSORT => q{ tmpfile=%s/bamsort_%s inputformat=sam verbose=0 index=1 md5=1 md5filename=%s.md5 indexfilename=%s.bai O=%s};
# out_bamname, out_bamname, out_bamname

#genome.fa.fai gc_windows.bed[.gz] in.bam out_path [chr_idx]
const my $BRASS_COVERAGE => q{ %s %s %s %s};

const my $NORM_CN_R => q{normalise_cn_by_gc_and_fb_reads.R %s %s %s %s %s %s %s};
const my $MET_HAST_R => q{metropolis_hastings_inversions.R %s %s %s};
const my $DEL_FB_R => q{filter_small_deletions_and_fb_artefacts.R %s %s %s %s};

const my $ISIZE_CHR => 5;
const my $CLEAN_ISIZE => q{%s view -f 33 -F 3868 %s %s | %s %s > %s};

## group
const my $BRASS_GROUP => q{ -I %s};
# extreme depth, repeats.gz, tumour.brm.bam, normal.brm.bam[, normal_panel.bam]
const my $FILTER_GROUP => q{ -i - -t %s -o %s};
# tumour name
const my $REDIR_GROUP => q{(%s) > %s};
#/software/CGP/projects/brass/bin/brass-group  | /nfs/users/nfs_k/kr2/git/brass/perl/bin/filter_groups.pl -t PD3904a > output.rearr

## filter
const my $BRASS_FILTER => q{ -seq_depth 25.1 -blat %s -ref %s -tumour %s -infile %s -outfile %s -min_tumour_count_high %d};
# path/to/blat, genome.fa, tumour name, groups_in, groups_out
#perl ~kr2/git/brass/perl/bin/brassI_filter.pl

## assemble
const my $BRASS_ASSEMBLE => q{ -X -m mem -O bedpe -r %s -T %s -o %s %s %s %s};
# extreme depth, genome.fa, tmp, output.tab, groups, tumourbam, normalbam

## grass
const my $GRASS => q{ -genome_cache %s -ref %s -species %s -assembly %s -platform %s -protocol %s -tumour %s -normal %s -file %s -add_header 'brassVersion=%s'};

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

    my $rg_prefix;
    if($index == 1) {
      $rg_prefix = 'T';
    }
    else {
      $rg_prefix = 'N';
    }

    my $command = _which('samtools');
    $command .= ' view -F 3854 -q 6 -u '.$input;
    $command .= ' | ';
    $command .= _which('bedtools');
    $command .= sprintf $BED_INTERSECT_V, 'stdin', $options->{'depth'};
    $command .= ' | ';
    $command .= _which('bamcollate2');
    $command .= sprintf $BAMCOLLATE, $tmp, $index;
    $command .= ' | ';
    $command .= "$^X ";
    $command .= _which('brassI_prep_bam.pl');
    $command .= sprintf $PREP_BAM, $input.'.bas', $rg_prefix;
    $command .= " -i ".join(',', valid_seqs($options));
    $command .= ' | ';
    $command .= _which('bamsort');
    $command .= sprintf $BAMSORT, $tmp, $index, $outbam, $outbam, $outbam;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), qq{bash -c 'set -o pipefail; $command'}, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub cover {
 my ($index_in, $options) = @_;

  my $tmp = $options->{'tmp'};

	# first handle the easy bit, skip if limit not set
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

  my @seqs = Sanger::CGP::Brass::Implement::valid_seqs($options);
  my @indicies = limited_indices($options, $index_in, $options->{'fai_count'});

  my $tmp_cover = File::Spec->catdir($tmp, 'cover');
  make_path($tmp_cover) unless(-e $tmp_cover);

  for my $index(@indicies) {
    # warning, extra layer of looping here, but individual success files still
    for my $samp_type(qw(tumour normal)) {
      next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $samp_type, $index);

      my $command = "$^X ";
      $command .= _which('compute_coverage.pl');
      $command .= sprintf $BRASS_COVERAGE, $options->{'gcbins'}, $options->{$samp_type}, $tmp_cover, $seqs[$index-1];

      PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $samp_type, $index);

      PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $samp_type, $index);
    }
  }
}

sub merge {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tmp_cover = File::Spec->catdir($tmp, 'cover');

  my $command_t = "$^X ";
  $command_t .= _which('coverage_merge.pl');
  $command_t .= sprintf ' %s %s %s', $options->{'genome'}.'.fai', $options->{'safe_tumour_name'}, $tmp_cover;
  $command_t .= ' '.join(',', valid_seqs($options));

  my $command_n = "$^X ";
  $command_n .= _which('coverage_merge.pl');
  $command_n .= sprintf ' %s %s %s', $options->{'genome'}.'.fai', $options->{'safe_normal_name'}, $tmp_cover;
  $command_n .= ' '.join(',', valid_seqs($options));

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), [$command_t, $command_n], 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub normcn {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tmp_cover = File::Spec->catdir($tmp, 'cover');
  my $tum_fb = sprintf '%s/%s.ngscn.fb_reads.bed.gz', $tmp_cover, $options->{'safe_tumour_name'};
  my $norm_fb = sprintf '%s/%s.ngscn.fb_reads.bed.gz', $tmp_cover, $options->{'safe_normal_name'};
  my $normcn_stub = sprintf '%s/%s_vs_%s.ngscn', $tmp_cover, $options->{'safe_tumour_name'}, $options->{'safe_normal_name'};

  my $command = _which('Rscript').' ';
  $command .= _Rpath().'/';
  $command .= sprintf $NORM_CN_R, $tum_fb,
                                  $norm_fb,
                                  $options->{'Ploidy'},
                                  $options->{'Acf'},
                                  $normcn_stub,
                                  $options->{'centtel'},
                                  $options->{'cytoband'};

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  check_tsv_file($normcn_stub.'.abs_cn.bg', 5);
  check_tsv_file($normcn_stub.'.segments.abs_cn.bg', 6);


  my $r_stub = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'});

  move($normcn_stub.'.abs_cn.bg',$r_stub.'.ngscn.abs_cn.bg') || die "Move failed: $!\n";
  move($normcn_stub.'.segments.abs_cn.bg', $r_stub.'.ngscn.segments.abs_cn.bg') || die "Move failed: $!\n";
  if(-e $normcn_stub.'.diagnostic_plots.pdf') {
    move($normcn_stub.'.diagnostic_plots.pdf', $r_stub.'.ngscn.diagnostic_plots.pdf') || die "Move failed: $!\n";
  }

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub check_tsv_file {
  my ($file, $expect_cols) = @_;
  open my $CHK, '<', $file || die $!;
  while(<$CHK>) {
    next if($_ =~ m/^#/);
    my @F = split /\t/, $_;
    die "FATAL: Row $. of $file does not have $expect_cols\n" if(scalar @F != $expect_cols);
  }
  close $CHK;
}

sub group {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups');
  my $tumour = File::Spec->catfile($tmp, $options->{'safe_tumour_name'}).'.brm.bam';
  my $normal = File::Spec->catfile($tmp, $options->{'safe_normal_name'}).'.brm.bam';

  my $command = _which('brass-group');
  $command .= sprintf $BRASS_GROUP, $options->{'depth'};
  $command .= ' -n '.$options->{'mingroup'};
  $command .= ' -F '.$options->{'repeats'} if(exists $options->{'repeats'});
  $command .= " $tumour $normal";
  $command .= ' | ';
  $command .= "$^X ";
  $command .= _which('brassI_pre_filter.pl');
  $command .= sprintf $FILTER_GROUP, $options->{'tumour_name'}, $groups;
  $command .= " -n $options->{filter}" if(exists $options->{'filter'});

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), qq{bash -c 'set -o pipefail; $command'}, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub isize {
  my $options = shift;
  my $tmp = $options->{'tmp'};

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
  my($chr_isize)=_getChr("$decomp.fai");

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tumour_isize = $options->{'outdir'}.'/'.$options->{'safe_tumour_name'}.'.insert_size_distr';

  my $command = sprintf $CLEAN_ISIZE, _which('samtools'), $options->{'tumour'}, $chr_isize, $^X, _which('corrected_insertsize.pl'), $tumour_isize;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), qq{bash -c 'set -o pipefail; $command'}, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}


sub _getChr {
	my($fai)=shift;

  my $fai_seqs = capture_stdout { system('cut', '-f', 1, $fai ); };
  my @array = split /\n/, $fai_seqs;
  if( grep( /^$ISIZE_CHR$/, @array ) ) {
     return $ISIZE_CHR;
	}elsif( grep(/^chr$ISIZE_CHR$/, @array ) ) {
		  return "chr$ISIZE_CHR";
  }else{
		return $array[4];
	}
}




sub filter {
  my $options = shift;

  my $tmp = $options->{'tmp'};

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
  my $brassI_filter = "$^X ";
  $brassI_filter .= _which('brassI_filter.pl');
  $brassI_filter .= sprintf $BRASS_FILTER, $blat,
                                      $decomp,
                                      $options->{'safe_tumour_name'},
                                      $groups,
                                      $filtered,
                                      $options->{'minkeep'};
  $brassI_filter .= ' -augment_off';

  my $bedpe = $filtered.'.bedpe';

  my $bedpe_no_head = $bedpe.'nohead';
  my $rem_head = "grep -v '^#' $bedpe > $bedpe_no_head";

  my $r_stub = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'});

  my $met_hasting = _which('Rscript').' ';
  $met_hasting .= _Rpath().'/';
  $met_hasting .= sprintf $MET_HAST_R, $bedpe_no_head, $options->{'safe_tumour_name'}, $r_stub;

  my $isize_dist = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'.insert_size_distr');
  my $fb_art = $r_stub.'.is_fb_artefact.txt';

  my $rfilt = $r_stub.'.r2';

  my $del_fb = _which('Rscript').' ';
  $del_fb .= _Rpath().'/';
  $del_fb .= sprintf $DEL_FB_R, $bedpe_no_head, $isize_dist, $fb_art, $rfilt;

  my $merge_file = $r_stub.'.r3';
  my $merge_grps = $^X.' '._which('merge_double_rgs.pl');
  $merge_grps .= " $rfilt $merge_file";

  my $abs_bkp_file = $r_stub.'.r4';
  my $abs_bkp = $^X.' '._which('get_abs_bkpts_from_clipped_reads.pl');
  $abs_bkp .= ' -fasta '.$options->{'genome'};
  $abs_bkp .= ' -out '.$abs_bkp_file;
  $abs_bkp .= ' '.$options->{'tumour'};
  $abs_bkp .= ' '.$merge_file;


  my $score_file = $r_stub.'.r5.scores';
  my $remap_file = $r_stub.'.r5';
  my $tumour_brm = File::Spec->catfile($tmp, sanitised_sample_from_bam($options->{'tumour'})).'.brm.bam';
  my $remap_micro = $^X.' '._which('filter_with_microbes_and_remapping.pl');
  $remap_micro .= sprintf ' -virus_db %s -bacterial_db_stub %s -scores_output_file %s -tmpdir %s -score_alg %s -search_cores %d -groups_file %s',
                          $options->{'viral'},
                          $options->{'microbe'},
                          $score_file,
                          File::Spec->catdir($tmp,'remap_micro'),
                          'ssearch36',
                          $options->{'threads'},
                          $bedpe_no_head;
  $remap_micro .= sprintf ' %s %s %s %s',
                          $abs_bkp_file,
                          $tumour_brm,
                          $options->{'genome'},
                          $remap_file;

  my $abs_cn = $r_stub.'.ngscn.abs_cn.bg';
  my $seg_cn = $r_stub.'.ngscn.segments.abs_cn.bg';

  my $rg_cns = _which('Rscript').' ';
  $rg_cns .= _Rpath().'/get_rg_cns.R';
  $rg_cns .= ' '.$remap_file;
  $rg_cns .= ' '.$abs_cn;
  $rg_cns .= ' '.$seg_cn;
  $rg_cns .= ' '.$options->{'tumour'};
  $rg_cns .= ' '.$options->{'Acf'};
  $rg_cns .= ' '.$options->{'centtel'};
  $rg_cns .= ' '.$options->{'GenderChr'};
  $rg_cns .= ' '.$options->{'GenderChrFound'};
  $rg_cns .= ' '.$tmp;

  my $match_lib_file = $r_stub.'.r6';
  my $filtered_cn = $r_stub.'.cn_filtered';
  my $match_lib = $^X.' '._which('match_rg_patterns_to_library.pl');
  $match_lib .= ' -filtered_bedpe '.$match_lib_file;
  $match_lib .= ' -acf '.$options->{'Acf'};
  $match_lib .= ' -ploidy '.$options->{'Ploidy'};
  $match_lib .= ' -min_cn_change '.$options->{'mincn'};
  $match_lib .= ' -filt_cn_out '.$filtered_cn;
  $match_lib .= " $remap_file";
  $match_lib .= ' '.$r_stub.'.ngscn.abs_cn.bg.rg_cns';

  my $final_file = $r_stub.'.groups.clean.bedpe';
  my $final = $^X.' '._which('collate_rg_regions.pl');
  $final .= ' '.$r_stub.'.groups.filtered.bedpe';
  $final .= ' '.$remap_file;
  $final .= ' '.$match_lib_file;
  $final .= ' '.$final_file;

  my $add_blat = "$^X ";
  $add_blat .= _which('brassI_filter.pl');
  $add_blat .= sprintf $BRASS_FILTER, $blat,
                                      $decomp,
                                      $options->{'safe_tumour_name'},
                                      $final_file,
                                      $final_file,
                                      $options->{'minkeep'};
  $add_blat .= ' -augment_only';

  #[$brassI_filter, $met_hasting, $del_fb, $merge_grps, $abs_bkp, $remap_micro, $rg_cns, $match_lib, $final]
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'brassI_filter')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $brassI_filter, 'brassI_filter');
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'brassI_filter');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'rem_head')) {
		my ($sto, $ste, $exit) = capture{ system([0,1],$rem_head) };
  	die "ERROR executing: $rem_head\n" if($exit > 1);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'rem_head');
    sleep 5;
  }

  my $group_count = count_lines($bedpe_no_head);

  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'met_hasting')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $met_hasting, 'met_hasting') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'met_hasting');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'del_fb')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $del_fb, 'del_fb') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'del_fb');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_grps')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $merge_grps, 'merge_grps') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_grps');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'abs_bkp')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $abs_bkp, 'abs_bkp') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'abs_bkp');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'remap_micro')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $remap_micro, 'remap_micro') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'remap_micro');
    sleep 5;
  }

  my $r4_lines = count_lines($abs_bkp_file, 1);
  my $r5scr_lines = count_lines($score_file, 1);

  die "Line count mismatch: $abs_bkp_file ($r4_lines) vs $score_file ($r5scr_lines)\n" if($r4_lines != $r5scr_lines);

  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'rg_cns')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $rg_cns, 'rg_cns') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'rg_cns');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'match_lib')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $match_lib, 'match_lib') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'match_lib');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'final')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $final, 'final') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'final');
    sleep 5;
  }
  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'add_blat')) {
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $add_blat, 'add_blat') if($group_count);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'add_blat');
    sleep 5;
  }
  if($group_count == 0) {
  	copy($bedpe, $final_file) or die "Failed to copy: $!";
  }
}

sub count_lines {
  my ($file, $allow_no_file) = @_;
  return 0 if($allow_no_file && ! -e $file);
  open my $FH, '<', $file or die $!;
  while(<$FH>){}
  my $lines = $.;
  close $FH;
  return $lines;
}

sub split_filtered {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups.clean.bedpe');

  my $split_dir = File::Spec->catdir($tmp, 'split');
  remove_tree($split_dir) if(-d $split_dir);
  make_path($split_dir);

  # before running this we need to know if there are any non '#' lines
  my ($sto, $ste, $exit) = capture { system([0,1], "grep -vc '^#' $groups"); };

	if($exit == 0 && $sto > 0) {
		my $command = sprintf q{grep -v "^#" %s | split --suffix-length=7 --numeric-suffixes --verbose --lines=%s - %s/split.},
													$groups, $ASSEMBLE_SPLIT, $split_dir;

		PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), qq{bash -c 'set -o pipefail; $command'}, 0);
  }
  elsif($exit == 1 && $sto == 0) {
	  my $splitOne = $split_dir.'/split.0000000';
    `touch $splitOne`;
    warn "No data survived filtering (split_filtered)\n";
  }
  else {
    die "ERROR: Unexpected result when parsing non-header lines from $groups\n";
  }

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
    $split_file .= sprintf '%07d', $index-1;

    # changed to be a real tmpdir as we always run with -X
    my $tmp_assemble = tempdir( 'brassAssembly_XXXXXX', TMPDIR => 1, CLEANUP => 1 );

    my $assemble_out = File::Spec->catdir($tmp, 'assemble');
    make_path($assemble_out) unless(-e $assemble_out);


    my $assembled = File::Spec->catfile($assemble_out, 'bedpe.');
    $assembled .= sprintf '%07d', $index-1;

    my $command = "$^X ";
    $command .= _which('brass-assemble');
    $command .= sprintf $BRASS_ASSEMBLE, $options->{'genome'},
                                          $tmp_assemble,
                                          $assembled,
                                          $split_file,
                                          $options->{'tumour'},
                                          $options->{'normal'};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub grass {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $assembled = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.assembled.bedpe');
  my $groups = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.groups.clean.bedpe');
  my $merge = sprintf '(cat %s/assemble/bedpe.* | sort -k1,1 -k 2,2n > %s)', $tmp, $assembled;

  my $tumour = File::Spec->catfile($tmp, $options->{'safe_tumour_name'}).'.brm.bam';
  my $normal = File::Spec->catfile($tmp, $options->{'safe_normal_name'}).'.brm.bam';

  my $annot_phaseI_prefix = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'_ann.groups.clean');
  my $annot_phaseII_prefix = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'_ann.assembled');

  my $final = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.annot');

  my $combine_cmd = "$^X ";
  $combine_cmd .= _which('combineResults.pl');
  $combine_cmd .= ' '.$annot_phaseI_prefix;
  $combine_cmd .= ' '.$annot_phaseII_prefix;
  $combine_cmd .= ' '.$final;
  $combine_cmd .= ' '.$options->{'Ploidy'};
  $combine_cmd .= ' '.$options->{'Acf'};
  $combine_cmd .= ' '.$options->{'PloidyAcfState'};

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'),
                                            [ qq{bash -c 'set -o pipefail; $merge'},
                                              _grass($options, $assembled),
                                              _grass($options, $groups),
                                              $combine_cmd],
                                            0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return 1;
}

sub _grass {
  my ($options, $input) = @_;
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
                              $input,
                              $VERSION;
  return $grass_cmd;
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
  my ($options, $type) = @_;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $type);

  my $annotated = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.annot.'.$type);
  my $sorted = $annotated.'.srt';

  my $tabix_type = 'vcf';
  my $sort_str = '-k 1,1 -k 2,2n';
  if($type eq 'bedpe') {
    $tabix_type = 'bed';
    $sort_str = '-k 1,1 -k 2,3n';
  }

  my $header = sprintf q{(grep '^#' %s > %s)}, $annotated, $sorted;
  my $sort = sprintf q{(grep -v "^#" %s | sort %s >> %s)}, $annotated, $sort_str, $sorted;

  my $type_gz = $annotated.'.gz';
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $sorted, $type_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p %s %s', $tabix_type, $type_gz;

  # before running this we need to know if there are any non '#' lines
  my ($sto, $ste, $exit) = capture { system([0,1], "grep -vc '^#' $annotated"); };

  my @commands = ($header);

  if($exit == 0 && $sto > 0) {
    push @commands, qq{bash -c 'set -o pipefail; $sort'};
  }
  elsif($exit == 1 && $sto == 0) {
    warn "No data survived (tabix, $type)\n";
  }
  else {
    die "ERROR: Unexpected result when parsing non-header lines from $annotated\n";
  }

  push @commands, $bgzip, $tabix;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, $type);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $type);
  return 1;
}

sub bedGraphToBigWig {
  my ($options, $bedgraph) = @_;
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  # required file
  my $chr_sizes = File::Spec->catfile($tmp, 'chrSizes.txt');
  my $chrCmd = sprintf 'cut -f 1,2 %s > %s', $options->{'genome'}.'.fai', $chr_sizes;

  # Need to strip down the cn file as it's not really a bedgraph
  my $realBg = File::Spec->catfile($tmp, 'real.bg');
  my $bgCmd = sprintf 'cut -f 1-4 %s > %s', $bedgraph, $realBg;

  my $bw = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'}.'.ngscn.abs_cn.bw');
  my $command = sprintf '%s -i %s -c %s -o %s', _which('bg2bw'), $realBg, $chr_sizes, $bw;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), [$chrCmd, $bgCmd, $command], 0);

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

sub _Rpath {
  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Brass::Implement') unless(-e File::Spec->catdir($mod_path, 'Rscripts'));

  my $rpath = File::Spec->catdir($mod_path, 'Rscripts');
  return $rpath;
}

sub get_ascat_summary {
  my ($options) = @_;
  my $ret_val = 0;
  $options->{'PloidyAcfState'} = q{};
  if(defined $options->{'ascat_summary'}) {
    $ret_val = 1;
    open my $SUMM, '<', $options->{'ascat_summary'};
    while(my $line = <$SUMM>) {
      chomp $line;
      my ($key, $value) = split /[[:space:]]+/, $line;
      next unless($key eq 'rho' || $key eq 'Ploidy' || $key eq 'GenderChr' || $key eq 'GenderChrFound');
      $key = 'Acf' if($key eq 'rho');
      $options->{$key} = $value;
    }

    if($options->{'Ploidy'} eq '?') {
      $options->{'PloidyAcfState'} = 'FAIL-DEFAULT'
    }
    else {
      $options->{'PloidyAcfState'} = 'REAL'
    }
  }

  if($options->{'PloidyAcfState'} ne 'REAL') {
    $options->{'Ploidy'} = 2;
    $options->{'Acf'} = 0.75;
    $options->{'PloidyAcfState'} = 'DEFAULT';
    $options->{'GenderChr'} = 'Y' unless(defined $options->{'GenderChr'});
    $options->{'GenderChrFound'} = 'Y' unless(defined $options->{'GenderChrFound'});
  }
  return $ret_val;
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_\-.]/_/ig; # sanitise sample name
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
      $tum_plat = $_->{'PL'};
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
      $norm_plat = $_->{'PL'};
    }
  }

  my $platform;
  if(defined $tum_plat && defined $norm_plat) {
   	$tum_plat = uc $tum_plat;
	  $norm_plat = uc $norm_plat;
    die "BAMs have different sequencing platform: $tum_plat, $norm_plat\n" if($norm_plat ne $tum_plat);
    $platform = $norm_plat;
  }

  die "Specified tumour name doesn't match BAM header" if(defined $options->{'tumour_name'} && defined $tum_name && $options->{'tumour_name'} ne $tum_name);
  die "Specified normal name doesn't match BAM header" if(defined $options->{'normal_name'} && defined $norm_name && $options->{'normal_name'} ne $norm_name);

  $options->{'tumour_name'} = $tum_name unless(defined $options->{'tumour_name'});
  $options->{'normal_name'} = $norm_name unless(defined $options->{'normal_name'});

  if(defined $options->{'platform'} && defined $platform && $options->{'platform'} ne $platform) {
  	warn "Specified platform name doesn't match BAM headers, default to specifed";
  }
  else {
  	$options->{'platform'} = $platform unless(defined $options->{'platform'});
  }

  return $options;
}

sub valid_seqs {
  my $options = shift;
  my @cent_seqs;
  open my $CENT, '<', $options->{'centtel'} or die $!;
  while (my $line = <$CENT>) {
    my ($item) = $line =~ m/^([^\t]+)\t/;
    next if($item eq 'chr');
    push @cent_seqs, $item;
  }
  close $CENT;

  my @good_seqs;
  open my $FAI, '<', $options->{'genome'}.'.fai' or die $!;
  while(my $line = <$FAI>) {
    my $chr = (split /\t/, $line)[0];
    push @good_seqs, $chr if(first {$chr eq $_} @cent_seqs);
  }
  close $FAI;

  return @good_seqs;
}

1;
