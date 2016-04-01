#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014,2015 Genome Research Ltd.
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


BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage qw(pod2usage);

use Bio::Brass::Merge;

use PCAP::Cli;

my $options = &setup;

my $merger = Bio::Brass::Merge->new(comment_char => '#',
                                    normal_groups => $options->{'normals'},
                                    analysis_groups => $options->{'input'},
                                    tumour => $options->{'tumour'},
                                    );

my $ofh;
if($options->{'output'} eq '-') {
  $ofh = *STDOUT;
}
else {
  open $ofh, '>', $options->{'output'};
}

print $ofh $merger->merge_headers or die "Failed to write to: $options->{output}";

$merger->merge_records($ofh);

close $ofh if($options->{'output'} ne '-');






sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'i|input=s' => \$opts{'input'},
              'o|output=s' => \$opts{'output'},
              'n|normals=s' => \$opts{'normals'},
              't|tumour=s' => \$opts{'tumour'},
              'v|version' => \$opts{'version'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ', Bio::Brass->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::file_for_reading('input', $opts{'input'}) unless($opts{'input'} eq q{-});
  PCAP::Cli::file_for_reading('normals', $opts{'normals'}) if(defined $opts{'normals'});

  return \%opts;
}

__END__

=head1 brassI_pre_filter.pl

Filters the raw groups against groups file generated from a normal panel and other low-level cleanup.

=head1 SYNOPSIS

brassI_pre_filter.pl [options]

  Required parameters:
    -output    -o   File or '-' for STDOUT
    -input     -i   Brass groups file or '-' for SDTIN
    -normals   -n   bgzip tabix-ed normal panel groups file
    -tumour    -t   Tumour sample name

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.
    brassI_pre_filter.pl -i Tumour_vs_Normal.groups -o data.filtered.groups -n brass_np.groups.gz -t Tumour
