#!/usr/bin/perl

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


use FindBin;
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
const my @VALID_PROCESS => qw(input group filter split assemble grass tabix);
my %index_max = ( 'input'   => 2,
                  'group'  => 1,
                  'filter' => 1,
                  'split' => 1,
                  'assemble'   => -1,
                  'grass' => 1,
                  'tabix' => 1);

{
  my $options = setup();
  Sanger::CGP::Brass::Implement::prepare($options);
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register any process that can run in parallel here
  $threads->add_function('input', \&Sanger::CGP::Brass::Implement::input);
  $threads->add_function('assemble', \&Sanger::CGP::Brass::Implement::assemble);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  $threads->run(2, 'input', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'input');

  Sanger::CGP::Brass::Implement::group($options) if(!exists $options->{'process'} || $options->{'process'} eq 'group');

  Sanger::CGP::Brass::Implement::filter($options) if(!exists $options->{'process'} || $options->{'process'} eq 'filter');

  Sanger::CGP::Brass::Implement::split_filtered($options) if(!exists $options->{'process'} || $options->{'process'} eq 'split');

  if(!exists $options->{'process'} || $options->{'process'} eq 'assemble') {
    $options->{'splits'} = Sanger::CGP::Brass::Implement::split_count($options);
    $threads->run($options->{'splits'}, 'assemble', $options)
  }

  Sanger::CGP::Brass::Implement::grass($options) if(!exists $options->{'process'} || $options->{'process'} eq 'grass');

  if(!exists $options->{'process'} || $options->{'process'} eq 'tabix') {
    Sanger::CGP::Brass::Implement::tabix($options);
    cleanup($options);
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  my $basefile = File::Spec->catfile($options->{'outdir'}, $options->{'safe_tumour_name'}.'_vs_'.$options->{'safe_normal_name'});

  unlink "$basefile.groups";
  unlink "$basefile.assembled.bedpe";
  unlink "$basefile.annot.vcf";
  unlink "$basefile.annot.vcf.srt";

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
              'a|ascat=s' => \$opts{'ascat'},
              'f|filter=s' => \$opts{'filter'},
              'e|exclude=s' => \$opts{'exclude'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'v|version' => \$opts{'version'},
              's|species=s' => \$opts{'species'},
              'as|assembly=s' => \$opts{'assembly'},
              'pr|protocol=s' => \$opts{'protocol'},
              'tn|tum_name=s' => \$opts{'tumour_name'},
              'nn|norm_name=s' => \$opts{'normal_name'},
              'pl|platform=s' => \$opts{'platform'},
              'gc|g_cache=s' => \$opts{'g_cache'},
              'l|limit=i' => \$opts{'limit'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ',Sanger::CGP::Pindel::Implement->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('depth', $opts{'depth'});
  PCAP::Cli::file_for_reading('genome', $opts{'genome'});
  PCAP::Cli::file_for_reading('repeats', $opts{'repeats'});
  PCAP::Cli::file_for_reading('g_cache', $opts{'g_cache'});
  PCAP::Cli::file_for_reading('ascat', $opts{'ascat'}) if(defined $opts{'ascat'});
  PCAP::Cli::file_for_reading('filter', $opts{'filter'}) if(defined $opts{'filter'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'ascat'} unless(defined $opts{'ascat'});
  delete $opts{'exclude'} unless(defined $opts{'exclude'});
  delete $opts{'filter'} unless(defined $opts{'filter'});
  delete $opts{'limit'} unless(defined $opts{'limit'});

  die "ERROR: No '.bas' file (with content) for tumour bam has been found $opts{tumour}\n" unless(-e $opts{'tumour'}.'.bas' && -s _ > 0);
  die "ERROR: No '.bas' file (with content) for normal bam has been found $opts{normal}\n" unless(-e $opts{'normal'}.'.bas' && -s _ > 0);

  for(@REQUIRED_PARAMS) {
    pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
  }

  # now safe to apply defaults
  $opts{'threads'} = 1 unless(defined $opts{'threads'});

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
        else {
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

  Sanger::CGP::Brass::Implement::get_bam_info(\%opts);
  die "Unable to find norm_name in BAM header please specify as option" unless(defined $opts{'normal_name'});
  die "Unable to find tum_name in BAM header please specify as option" unless(defined $opts{'tumour_name'});
  die "Unable to find platform in BAM header please specify as option" unless(defined $opts{'platform'});

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
    -repeats   -r   Repeat file, see 'make-repeat-file'
    -genome    -g   Genome fasta file
    -species   -s   Species name
    -assembly  -as  Assembly name
    -protocol  -pr  Sequencing protocol (WGS|WXS|RNA)
    -g_cache   -gc  Genome cache file.

  Optional
    -ascat     -a   ASCAT copynumber summary
    -platform  -pl  Sequencing platform (when not defined in BAM header)
    -tum_name  -tn  Tumour sample name (when not defined in BAM header)
    -norm_name -nn  Normal sample name (when not defined in BAM header)
    -exclude   -e   Exclude this list of ref sequences from processing, wildcard '%'
    -filter    -f   bgzip tabix-ed normal panel groups file
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 2 during 'input' process.

  Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.
    brass.pl -c 2 -o myout -t tumour.bam -n normal.bam

  Run with '-m' for possible input file types.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  input
  group
  filter
  assemble
  grass
  tabix

=item B<-index>

Possible index ranges for processes above are:

  input     = 1..2
  group     = 1
  filter    = 1
  split     = 1
  assemble  = ?
  grass     = 1
  tabix     = 1

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=back
