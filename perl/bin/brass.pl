#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2020 Genome Research Ltd.
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


use FindBin qw($Bin);
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy;

use PCAP::Cli;
use Sanger::CGP::Brass::Implement;

const my @REQUIRED_PARAMS => qw(species assembly protocol);
const my @VALID_PROCESS => qw(input cover merge normcn group isize filter split assemble grass tabix);
my %index_max = ( 'input'   => 2, # input and cover can run at same time
                  'cover' => -1,
                  'merge' => 1,
                  'group'  => 1,
                  'isize' => 1,
                  'normcn' => 1, # normcn, group and isize can run at same time
                  'filter' => 1,
                  'split' => 1,
                  'assemble'   => -1,
                  'grass' => 1,
                  'tabix' => 1);

{
  my $options = setup();
  my $threads = PCAP::Threaded->new($options->{'threads'});

  # register any process that can run in parallel here
  $threads->add_function('input', \&Sanger::CGP::Brass::Implement::input);
  $threads->add_function('cover', \&Sanger::CGP::Brass::Implement::cover);
  $threads->add_function('assemble', \&Sanger::CGP::Brass::Implement::assemble);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  $threads->run(2, 'input', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'input');

  if(!exists $options->{'process'} || $options->{'process'} eq 'cover') {
    $options->{'fai_count'} = scalar Sanger::CGP::Brass::Implement::valid_seqs($options);
    my $jobs = $options->{'fai_count'};
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'cover', $options);
  }

  Sanger::CGP::Brass::Implement::merge($options) if(!exists $options->{'process'} || $options->{'process'} eq 'merge');

  Sanger::CGP::Brass::Implement::group($options) if(!exists $options->{'process'} || $options->{'process'} eq 'group');

  Sanger::CGP::Brass::Implement::isize($options) if(!exists $options->{'process'} || $options->{'process'} eq 'isize');

  Sanger::CGP::Brass::Implement::normcn($options) if(!exists $options->{'process'} || $options->{'process'} eq 'normcn');

  Sanger::CGP::Brass::Implement::filter($options) if(!exists $options->{'process'} || $options->{'process'} eq 'filter');

  Sanger::CGP::Brass::Implement::split_filtered($options) if(!exists $options->{'process'} || $options->{'process'} eq 'split');

  if(!exists $options->{'process'} || $options->{'process'} eq 'assemble') {
    $options->{'splits'} = Sanger::CGP::Brass::Implement::split_count($options);
    my $jobs = $options->{'splits'};
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'assemble', $options);
  }

  Sanger::CGP::Brass::Implement::grass($options) if(!exists $options->{'process'} || $options->{'process'} eq 'grass');

  if(!exists $options->{'process'} || $options->{'process'} eq 'tabix') {
    Sanger::CGP::Brass::Implement::tabix($options, 'vcf');
    Sanger::CGP::Brass::Implement::tabix($options, 'bedpe');
    cleanup($options) unless($options->{'noclean'});
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  my $outdir = $options->{'outdir'};

  my $basefile = File::Spec->catfile($outdir, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'});

  # bigwig the copynumber file:
  Sanger::CGP::Brass::Implement::bedGraphToBigWig($options, $basefile.'.ngscn.abs_cn.bg');

  my $intdir = File::Spec->catdir($outdir, 'intermediates/.');
  make_path($intdir) unless(-d $intdir);

  my @to_gzip = ( $basefile.'.ngscn.abs_cn.bg',
  								$basefile.'.ngscn.abs_cn.bg.rg_cns',
                  $basefile.'.ngscn.segments.abs_cn.bg',
                  $basefile.'.groups',
                );

  for my $file(@to_gzip) {
    if(-e $file) {
      system(qq{gzip $file}) and die $!;
    }
    if(-e "$file.gz") {
      move("$file.gz", $intdir);
    }
  }

  # move the BRM files
  my @brm = glob qq{$tmpdir/*.brm.bam};
  for(@brm) {
    move($_, "$outdir/.");
  }

  # index the BRM files
  @brm = glob qq{$outdir/*.brm.bam};
  for(@brm) {
    system(qq{samtools index $_}) and die $!;
  }

  # read counts
  for my $sample($options->{'safe_tumour_name'}, $options->{'safe_normal_name'}) {
    for my $type(qw(.ngscn.bed.gz .ngscn.fb_reads.bed.gz)) {
      my $source = $tmpdir.'/cover/'.$sample.$type;
      move($source, $intdir) if(-e $source);
    }
  }

  my @base_delete = qw( _ann.assembled.bedpe
                        _ann.assembled.vcf
                        _ann.groups.clean.bedpe
                        _ann.groups.clean.vcf
                        .annot.vcf
                        .annot.vcf.srt
                        .annot.bedpe
                        .annot.bedpe.srt
                        .assembled.bedpe
                        .groups.filtered.bedpe.preclean
                        .groups.filtered.bedpenohead
                    );
  for my $to_del(@base_delete) {
    my $file = $basefile.$to_del;
    unlink $file if(-e $file);
  }
  my @base_move = glob qq{$basefile.r* $outdir/*.insert_size_distr $outdir/*.pdf $outdir/*.ascat*};
  push @base_move,  "$basefile.groups.filtered.bedpe",
                    "$basefile.is_fb_artefact.txt",
                    "$basefile.groups.clean.bedpe",
                    "$basefile.cn_filtered";

  for my $to_move(@base_move) {
    move($to_move, $intdir) if(-e $to_move);
  }

  # ascat files may not be in base output folder so copy from inputs if available (as can be optional)
  copy($options->{'ascat_summary'}, $intdir) if(-e $options->{'ascat_summary'});

	my $intermediates_dir = File::Spec->catdir($outdir, 'intermediates');
  my $targz = $basefile.'.intermediates.tar.gz';
  unless(-e $targz && -s $targz) {
  	system(qq{tar -C $outdir -zcf $targz intermediates}) and die $!;
  }

  my $tmplogs = File::Spec->catdir($tmpdir, 'logs');
  if(-e $tmplogs) {
    move($tmplogs, File::Spec->catdir($outdir, 'logs')) or die $!;
  }

	remove_tree $intermediates_dir if(-e $intermediates_dir);
  remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'c|cpus=i' => \$opts{'threads'},
              'o|outdir=s' => \$opts{'outdir'},
              't|tumour=s' => \$opts{'tumour'},
              'n|normal=s' => \$opts{'normal'},
              'd|depth=s' => \$opts{'depth'},
              'r|repeats=s' => \$opts{'repeats'},
              'g|genome=s' => \$opts{'genome'},
              'f|filter=s' => \$opts{'filter'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'b|gcbins=s' => \$opts{'gcbins'},
              'v|version' => \$opts{'version'},
              's|species=s' => \$opts{'species'},
              'ct|centtel=s' => \$opts{'centtel'},
              'cb|cytoband=s' => \$opts{'cytoband'},
              'ss|sampstat=s' => \$opts{'ascat_summary'},
              'as|assembly=s' => \$opts{'assembly'},
              'pr|protocol=s' => \$opts{'protocol'},
              'tn|tum_name=s' => \$opts{'tumour_name'},
              'nn|norm_name=s' => \$opts{'normal_name'},
              'pl|platform=s' => \$opts{'platform'},
              'gc|g_cache=s' => \$opts{'g_cache'},
              'vi|viral=s' => \$opts{'viral'},
              'mi|microbe=s' => \$opts{'microbe'},
              'cn|mincn=s' => \$opts{'mincn'},
              'j|mingroup=i' => \$opts{'mingroup'},
              'k|minkeep=i' => \$opts{'minkeep'},
##    -minkeep   -k   Minmum reads to retain a group [4]. ## disabled until metropolis_hastings_inversions.R can handle variable
              'l|limit=i' => \$opts{'limit'},
              'x|noclean' => \$opts{'noclean'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1, -exit=>0) if(defined $opts{'h'});
  pod2usage(-verbose => 2, -exit=>0) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ',Sanger::CGP::Brass::Implement->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
  $opts{'outdir'} = File::Spec->rel2abs( $opts{'outdir'} );
  $opts{'outdir'} = File::Spec->catdir(File::Spec->curdir(), $opts{'outdir'}) unless($opts{'outdir'} =~ m|^/|);

  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('depth', $opts{'depth'});
  PCAP::Cli::file_for_reading('genome', $opts{'genome'});
  PCAP::Cli::file_for_reading('viral', $opts{'viral'});
  PCAP::Cli::file_for_reading('centtel', $opts{'centtel'});
  PCAP::Cli::file_for_reading('repeats', $opts{'repeats'}) if(defined $opts{'repeats'});
  PCAP::Cli::file_for_reading('g_cache', $opts{'g_cache'});
  PCAP::Cli::file_for_reading('filter', $opts{'filter'}) if(defined $opts{'filter'});
  PCAP::Cli::file_for_reading('cytoband', $opts{'cytoband'});
  PCAP::Cli::file_for_reading('ascat_summary', $opts{'ascat_summary'});

  for my $item(qw(tumour normal depth genome viral repeats g_cache filter ascat_summary centtel )) {
    $opts{$item} = File::Spec->rel2abs( $opts{$item} ) if(defined $opts{$item});
  }

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'repeats'} unless(defined $opts{'repeats'});
  delete $opts{'filter'} unless(defined $opts{'filter'});
  delete $opts{'limit'} unless(defined $opts{'limit'});

  $opts{'mincn'} ||= 0.3;

  die "ERROR: No '.bas' file (with content) for tumour bam has been found $opts{tumour}\n" unless(-e $opts{'tumour'}.'.bas' && -s _ > 0);
  die "ERROR: No '.bas' file (with content) for normal bam has been found $opts{normal}\n" unless(-e $opts{'normal'}.'.bas' && -s _ > 0);

  die "ERROR: '-microbe' option has not been defined\n" unless(defined $opts{'microbe'});

  for(@REQUIRED_PARAMS) {
    pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
  }

  Sanger::CGP::Brass::Implement::get_bam_info(\%opts);
  die "Unable to find norm_name in BAM header please specify as option" unless(defined $opts{'normal_name'});
  die "Unable to find tum_name in BAM header please specify as option" unless(defined $opts{'tumour_name'});
  die "Unable to find platform in BAM header please specify as option" unless(defined $opts{'platform'});

  Sanger::CGP::Brass::Implement::prepare(\%opts);
  my $intermediates = File::Spec->catfile($opts{'outdir'}, $opts{'safe_tumour_name'}.'_vs_'.$opts{'safe_normal_name'}.'.intermediates.tar.gz');
  my $final_logs = File::Spec->catdir($opts{'outdir'}, 'logs');
  if(-e $intermediates && -e $final_logs) {
    warn "NOTE: Presence of intermediates.tar.gz and final logs directories suggests successful complete analysis, please delete to proceed: $intermediates\n";
    exit 0; # don't cause failure
  }

  # now safe to apply defaults
  $opts{'threads'} = 1 unless(defined $opts{'threads'});
  $opts{'mingroup'} = 2 unless(defined $opts{'mingroup'});
  $opts{'minkeep'} = 4 unless(defined $opts{'minkeep'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpBrass');
  make_path($tmpdir) unless(-d $tmpdir);
  $opts{'tmp'} = $tmpdir;

  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {

      my $max = $index_max{$opts{'process'}};
      if($max==-1){
        if(exists $opts{'limit'}) {
          $max = $opts{'limit'};
        }
        elsif($opts{'process'} eq 'cover') {
          $max = scalar Sanger::CGP::Brass::Implement::valid_seqs(\%opts);
        }
        elsif($opts{'process'} eq 'assemble') {
      	  $max = Sanger::CGP::Brass::Implement::split_count(\%opts);
      	}
      }
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      die "No max has been defined for this process type\n" if($max == 0);
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

	if(!defined $opts{'process'} || first { $opts{'process'} eq $_ } (qw(normcn filter grass))) {
  	if(Sanger::CGP::Brass::Implement::get_ascat_summary(\%opts)) {
    	die "Failed to find 'rho' in $opts{ascat_summary}\n" unless(exists $opts{'Acf'});
    	die "Failed to find 'Ploidy' in $opts{ascat_summary}\n" unless(exists $opts{'Ploidy'});
    	die "Failed to find 'GenderChr' in $opts{ascat_summary}\n" unless(exists $opts{'GenderChr'});
    	die "Failed to find 'GenderChrFound' in $opts{ascat_summary}\n" unless(exists $opts{'GenderChrFound'});
    }
  }

  return \%opts;
}


__END__

=head1 brass.pl

Reference implementation of Cancer Genome Project rearrangement calling pipeline.

=head1 SYNOPSIS

brass.pl [options]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -tumour    -t   Tumour BAM file
    -normal    -n   Normal BAM file
    -depth     -d   Regions of excessive sequencing depth to be ignored
    -genome    -g   Genome fasta file
    -species   -s   Species name
    -assembly  -as  Assembly name
    -protocol  -pr  Sequencing protocol (WGS only - remains for compat)
    -g_cache   -gc  Genome cache file.
    -viral     -vi  Virus sequences from NCBI
    -microbe   -mi  Microbe sequences file prefix from NCBI, please exclude '.N.fa.2bit'
    -gcbins    -b   5 column bed coord file, col 4 number of non-N bases, col 5 GC fraction.
    -cytoband    -cb  Cytoband file for a species build (can be obtained from UCSC)
    -centtel   -ct  TSV file of usable regions of the chromosome
                      Example in perl/share/Rdefault/
    -sampstat  -ss  ASCAT sample statistics file or file containing
                      rho 0.XXXXX [0.75] (~ 1-normalContamination)
                      Ploidy X.XXX [2.0]
                      GenderChr Y [Y]
                      GenderChrFound Y/N [Y]

  Optional
    -mingroup  -j   Minimum reads to call group [2].
    -mincn  	 -cn   Minimum CN change for copynumber_flag [0.3].
    -repeats   -r   Repeat file, see 'make-repeat-file' (legacy)
    -platform    -pl  Sequencing platform (when not defined in BAM header)

    -tum_name    -tn  Tumour sample name (when not defined in BAM header)
    -norm_name   -nn  Normal sample name (when not defined in BAM header)
    -filter      -f   bgzip tabix-ed normal panel groups file
    -noclean     -x   Don't tidyup the processing areas.
    -cpus        -c   Number of cores to use. [1]
                     - recommend max 2 during 'input' process.

  Targeted processing (further detail under PROCESS OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job
    -limit     -l   Define with -p and -i, see below

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.
    brass.pl -c 2 -o myout -t tumour.bam -n normal.bam

  Run with '-m' for possible input file types.

=head1 ONE SHOT PROCESSING

You can submit the whole analysis with a single command using mutliple cores on the same host.
Ensure that -p, -i and -l are not set when doing this.

=head1 PROCESS OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  input
  cover
  merge
  normcn
  group
  isize
  filter
  split
  assemble
  grass
  tabix

=item B<-index>

Possible index ranges for processes above are:

  input     = 1..2
  cover     = ? - use with -l
  merge     = 1
  normcn    = 1
  group     = 1
  isize     = 1
  filter    = 1
  split     = 1
  assemble  = ? - use with -l
  grass     = 1
  tabix     = 1

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=item B<-limit>

Allows -index of indeterminate size to executed on a known number of processes.

For example if 5 jobs are required to run cover step, setting limit to 2 will allow this to be run in 2 jobs:

 $ brass.pl ... -process cover -index 1 -limit 2
   # runs 1,3,5

 $ brass.pl ... -process cover -index 2 -limit 2
   # runs 2,4

=back
