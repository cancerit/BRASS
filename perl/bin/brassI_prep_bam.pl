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


# pre-filters brassI input to generate the correct headers in the bam file

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => 'all';

use Capture::Tiny qw { capture };
use PCAP::Bam::Bas;
use List::Util qw(first);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);

use Data::Dumper;

const my $MIN_MAPQ => 6;

my $options = &setup;
run($options);

sub run {
  my $options = shift;
  my @header;
  my $header_done = 0;
  my @good_seqs;
  my $p;
  $p = $options->{'prefix'} if(exists $options->{'prefix'});
  while(my $r1 = <>) {
    if($r1 =~ m/^\@/) {
      push @header, $r1;
      next;
    }
    unless($header_done) {
      my ($new_header, $all_seqs) = process_header($options, @header);
      @good_seqs = limit_seqs($options, $all_seqs);
      print $new_header;
      $header_done++;
      redo;
    }
    my $r2 = <>;
    my ($rname,$mapq,$rnext,$tlen,$seq) = (split /\t/, $r1)[2,4,6,8,9];
    next if($mapq < $MIN_MAPQ);
    next if($tlen != 0 && abs $tlen < (length $seq) * 2);
    next unless(first { $rname eq $_ } @good_seqs);
    next unless(first { $rnext eq $_ } @good_seqs);
    my $mapq2 = (split /\t/, $r2)[4];
    next if($mapq2 < $MIN_MAPQ);

    if($options->{'np'}) { # clear pbq to save space
      my @info = split /\t/, $r1;
      $info[10] = '*';
      $r1 = join qq{\t}, @info;
      @info = split /\t/, $r2;
      $info[10] = '*';
      $r2 = join qq{\t}, @info;
    }

    if($p) {
      $r1 =~ s/(\tRG:Z:)([^\t]+)/${1}${p}${2}/;
      $r2 =~ s/(\tRG:Z:)([^\t]+)/${1}${p}${2}/;
    }

    print $r1;
    print $r2;
  }
}

sub limit_seqs {
  my ($options, $raw_seqs) = @_;
  my @good_seqs = qw(= *);
  if(exists $options->{'exclude'}) {
    my @exclude = split /,/, $options->{'exclude'};
    my @exclude_patt;
    for my $ex(@exclude) {
      $ex =~ s/%/.+/;
      push @exclude_patt, $ex;
    }

    for my $sq(@{$raw_seqs}) {
      push @good_seqs, $sq unless(first { $sq =~ m/^$_$/ } @exclude_patt);
    }
  }
  else {
    push @good_seqs, @{$raw_seqs};
  }
  return @good_seqs;
}

#Once the new header is constructed (don't think this needs to be inside the script though):
#	samtools view -bF 3842 -q 1 ../PD13371a.bam | samtools reheader new_head.sam - > new.bam
#'-F' the flag values we don't want, see http://picard.sourceforge.net/explain-flags.html for what 3842 means
sub process_header {
  my ($options, @header) = @_;

  my @seq_names;
  my $new_header = '';

  # get readgroup ids and mean_insert_size out of bam.bas file
  my $bas_ob = PCAP::Bam::Bas->new($options->{'bas'});
  my $count_in_bas = get_bas_count($bas_ob);

  my $count_in_bam = 0;
  foreach my $row(@header) {
    if ($row !~ /^\@RG/) {
      $new_header .= $row;
      push @seq_names, $1 if($row =~ m/^\@SQ.*\tSN:([^\t]+)/);
      next;
    }

    # process the @RG line
    # order of lines in header MUST NOT CHANGE IN ANY WAY
    $count_in_bam++;
    my $rg_id = 0;
    $rg_id = $1 if($row =~ /\tID:([^\t]+)/);

    # die horribly if RG ids don't match
    unless ($bas_ob->get($rg_id, 'readgroup')) { die "Readgroup $rg_id in bam file but not in bas file $!"; }


    # The script needs use the information in the *.bas file to augment the Reagd-group line with the 'MI:' tag.
    #     This holds the 'maximum insert' size, this should be calculated as:
    #     mean_insert_size + (2* insert_size_sd)
    my $mean_insert_size = $bas_ob->get($rg_id, 'mean_insert_size');
    my $insert_size_sd = $bas_ob->get($rg_id, 'insert_size_sd');
    my $mi = int( $mean_insert_size + (2 * $insert_size_sd) );

    my $mi_field_found = 0;
    if($row =~ m/\tMI:(?:Z:)?[^\t]+/) {
      $row =~ s/\tMI:[^\t]+/\tMI:$mi/;
    }
    else {
      $row =~ s/\n/\tMI:$mi\n/;
    }

    if($row =~ m/\tPI:[^\t]*/) {
      $row =~ s/\tPI:[^\t]*/\tPI:$mean_insert_size/;
    }
    else {
      $row =~ s/\n/\tPI:$mean_insert_size\n/;
    }

    # prefix sample name with 'NP_'
    if($options->{'np'}) {
      my ($sample) = $row =~ m/\tSM:([^\t]+)/;
      $row =~ s/\tSM:[^\t]+/\tSM:NP_$sample/;
    }

    # now all correlation with BAS file complete we can modify the RGID if requested
    if(exists $options->{'prefix'}) {
      my $p = $options->{'prefix'};
      $row =~ s/(\tID:)([^\t]+)/${1}${p}${2}/;
    }

    $new_header .= $row;
  }

  # die horribly if number of RGs in bas and bam don't match
  unless ($count_in_bam == $count_in_bas) { die "Different number or readgroups in bam and bas files $!"; }

  return ($new_header, \@seq_names);
}

sub get_bas_count {
    my ($bas_ob) = @_;

    my @bas_rgs;
    my $method = 'read_groups';
    if ($bas_ob->can($method)) {
      @bas_rgs = $bas_ob->read_groups();
    }
    else {
      @bas_rgs = keys %{$bas_ob->{'_data'}};
    }
    return scalar @bas_rgs;
}

sub setup {
  my %opts;
  GetOptions( 'h|help'        => \$opts{'h'},
              'm|man'         => \$opts{'m'},
              'b|bas=s'       => \$opts{'bas'},
              'p|prefix:s'       => \$opts{'prefix'},
              'np|norm_panel' => \$opts{'np'},
              'e|exclude:s'   => \$opts{'exclude'},
  ) or pod2usage(1);

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  pod2usage(-message => qq{File not found: $opts{bas}}, -verbose => 1) unless(-e $opts{'bas'});
  pod2usage(-message => qq{Empty file: $opts{bas}}, -verbose => 1) unless(-s $opts{'bas'});

  delete $opts{'exclude'} unless(defined $opts{'exclude'});
  delete $opts{'prefix'} unless(defined $opts{'prefix'});

  return \%opts;
}


__END__

=head1 brassI_prep_bam.pl

Add max insert and exclude reads not mapped to expected refs.

=head1 SYNOPSIS

brassI_prep_bam.pl [options]

  Example
   ... | brassI_prep_bam.pl -b my.bam.bas -e NC_007605,hs37d5,GL% | <some digesting process>

=head1 OPTIONS

  Required parameters:
    -bas          -b    Bas statistics file for BAM being streamed

  Optional
    -prefix       -p    Prefix all readgroup IDs with this text to force unique between samples, (e.g. T, N)
    -exclude      -e    Exclude reads where self and mate are mapped to this ref name (or unmapped).
                         - csv, allows wild card of '%'
    -norm_panel   -np   For generation of normal panel input only

  Other:
    -help         -h    Brief help message.
    -man          -m    Full documentation.

=head1 Option Detail

=over 2

=item B<-bas>

Tab separated file of BAM file statistics generated by bam_stats.pl from PCAP-core.

Values for the following are required expected:

  mean_insert_size
  insert_size_sd

=item B<-exclude>

Ref names that should be excluded from prepared BAM.  Wildcard of '%' is allowed.

e.g.

  -e hs37d5,NC_007605,GL%

=item B<-norm_panel>

Use when preparing file for use in constructing a normal panel for filtering.

All sample names will be prefixed with 'NP_'.

All QUAL fields will be set to '*' to reduce output file size.

=back
