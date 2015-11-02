#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2015 Genome Research Ltd.
#
# Author: Yilong Li <yl3@sanger.ac.uk>
#
# This file is free software: you can redistribute it and/or modify it under
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

use FindBin;
use lib "$FindBin::Bin/../lib";

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use autodie qw(:all);
use Carp 'verbose';
use Scalar::Util qw(blessed);
use List::Util qw(sum);
$SIG { __DIE__ } = sub { Carp::confess(@_) };

use Sanger::CGP::CopyNumber::Segment;
use Sanger::CGP::CopyNumber::Segment::End;
use Sanger::CGP::CopyNumber::Segment::Genome;
use Sanger::CGP::Rearrangement::End;
use Sanger::CGP::Rearrangement::Path 'find_all_valid_paths';
use Sanger::CGP::Rearrangement::Group;
use Sanger::CGP::Rearrangement::ConnectedComponentSet;

#
# Handle options
#
pod2usage(-exitval => 1, -output => \*STDERR) if @ARGV < 2;

my %opts = (
    # Basic options
    acf    => 1,
    ploidy => 2,
    verbose => 0,
    help => 0,

    # Advanced options for extracting rearrangement patterns
    max_balanced_rg_dist => 500,
    max_balanced_rg_overlap => 100,
    max_foldback_distance => 5000,
    min_seg_size_for_cn => 5000,
    min_cn_change => 0.5,
    min_cn_bkpt_seg_size => 50000,
    keep_small_dels_and_tds => 1,
    max_depth => 3,
    max_shard_length => 10000,  ## TODO
    shard_bypassing_slop => 200,
    skip_chrs => '',

    within => 1e6,
    away_from => 3e6,
);
sub opt_handler {
    my ($opt_name, $opt_value) = @_;
    if ($opt_name eq "acf") {
        if ($opt_value < 0 || $opt_value > 100) {
            pod2usage("-acf value must be between 0 and 100!\n");
        }
        if ($opt_value > 1) { $opt_value /= 100; }
        $opts{acf} = $opt_value;
    }
    elsif ($opt_name eq "ploidy") {
        if ($opt_value <= 0) {
            pod2usage("-ploidy must be larger than 0!\n");
        }
        $opts{ploidy} = $opt_value;
    }
    elsif ($opt_name eq "max_balanced_rg_dist") {
        # if ($opt_value <= 0) {
        #     pod2usage("-max_balanced_rg_dist must be a non-negative integer!\n");
        # }
        $opts{max_balanced_rg_dist} = $opt_value;
    }
    elsif ($opt_name eq "max_foldback_distance") {
        if ($opt_value <= 0) {
            pod2usage("-max_foldback_distance must be a non-negative integer!\n");
        }
        $opts{max_foldback_distance} = $opt_value;
    }
    elsif ($opt_name eq "min_seg_size_for_cn") {
        if ($opt_value <= 0) {
            pod2usage("-min_seg_size_for_cn must be a non-negative integer!\n");
        }
        $opts{min_seg_size_for_cn} = $opt_value;
    }
    elsif ($opt_name eq "min_cn_change") {
        if ($opt_value <= 0) {
            pod2usage("-min_cn_change must be non-negative!\n");
        }
        $opts{min_cn_change} = $opt_value;
    }
    elsif ($opt_name eq "min_cn_bkpt_seg_size") {
        if ($opt_value <= 0) {
            pod2usage("-min_cn_bkpt_seg_size must be non-negative!\n");
        }
        $opts{min_cn_bkpt_seg_size} = $opt_value;
    }
    elsif ($opt_name eq "shard_bypassing_slop") {
        if ($opt_value <= 0) {
            pod2usage("-shard_bypassing_slop must be non-negative!\n");
        }
        $opts{shard_bypassing_slop} = $opt_value;
    }
    else {
        die;
    }
}
GetOptions(
    "acf=f" => \&opt_handler,
    "ploidy=f" => \&opt_handler,
    "verbose:1" => \$opts{verbose},
    "help" => \$opts{help},
    "max_balanced_rg_dist=i" => \&opt_handler,
    "max_foldback_distance=i" => \&opt_handler,
    "min_seg_size_for_cn=i" => \&opt_handler,
    "min_cn_change=f" => \&opt_handler,
    "min_cn_bkpt_seg_size=i" => \&opt_handler,
    "max_depth=i" => \$opts{max_depth},
    "filtered_bedpe=s" => \$opts{filtered_bedpe},
    "unique_bkpts=s" => \$opts{unique_bkpts},
    "naive_classification=s" => \$opts{naive_classification},
    "shard_bypassing_slop=i" => \&opt_handler,
    "skip_chrs=s" => \$opts{skip_chrs},
    "filt_cn_out=s" => \$opts{filt_cn_out},
);

my $rgs_file = shift;
my $cn_file  = shift;
pod2usage(-msg => "File '$rgs_file' doesn't exist!", -exitval => 2) if !-e $rgs_file;
pod2usage(-msg => "File '$cn_file' doesn't exist!", -exitval => 2) if !-e $cn_file;
pod2usage() if $opts{help};


#
# Some helper variables
#
my(@F);


#
# Actual commands for the code
#

# First read in all rearrangement and copy number data into a genome.
print STDERR "[" . `echo -n \`date\`` . "] Reading data in...\n" if $opts{verbose};
my $genome_of_cn_segs = Sanger::CGP::CopyNumber::Segment::Genome->new_from_BEDPE_and_CN_BED($rgs_file, $cn_file, %opts);

print STDERR "[" . `echo -n \`date\`` . "] Sorting segments on each chromosome...\n" if $opts{verbose};
$genome_of_cn_segs->sort_segments;

print STDERR "[" . `echo -n \`date\`` . "] Removing redundant breakpoints and filtering rearrangements...\n" if $opts{verbose};
$genome_of_cn_segs->preprocess_rearrangement_and_copy_number_data(%opts);

if ($opts{filtered_bedpe}) {
    print STDERR "[" . `echo -n \`date\`` . "] Outputting filtered rearrangements...\n" if $opts{verbose};
    open my $fh, '>', $opts{filtered_bedpe} || die $!;
    $genome_of_cn_segs->print_rearrangements_in_bedpe($fh, %opts);
    close $fh;
}
else {
    print STDERR "[" . `echo -n \`date\`` . "] Skipping rearrangement filtering outputting as -filtered_bedpe was not provided.\n" if $opts{verbose};
}

my $cn_fh = *STDOUT;
if(defined $opts{filt_cn_out}) {
  open $cn_fh, '>', $opts{filt_cn_out} || die $!;
}

$genome_of_cn_segs->print_rg_cns_bedpe($cn_fh);

close $cn_fh if(defined $opts{filt_cn_out});


#
# Subroutines
#
sub make_rg_groups_by_reciprocal_rg_recursion {
    ## Depth-first recursion. For each new component to be added, all the
    ## possible reciprocal rearrangement combinations are iterated through
    ## before passing to the next iteration of the recursion.

    my $remaining_components = shift;  ## These are remaining components,
                                       ## for which we have to find
                                       ## reciprocal component combinations.
    my $currently_included_components = shift;
    my $path = shift;
    die if Scalar::Util::blessed($remaining_components) ne 'Sanger::CGP::Rearrangement::ConnectedComponentSet';
    die if Scalar::Util::blessed($currently_included_components) ne 'Sanger::CGP::Rearrangement::ConnectedComponentSet';

    if ($remaining_components->size == 0) {
        my $new_rg_group = Sanger::CGP::Rearrangement::Group->new(
            components => [$currently_included_components->components_array],
            path => $path,
            target_component => $path->target_rg->component,
        );
        $new_rg_group->get_normalised_rg_patterns(%opts);  ## TODO
        return $new_rg_group;
    }
    else {
        my $new_remaining_components = Sanger::CGP::Rearrangement::ConnectedComponentSet->new(
            components => [$remaining_components->components_array]
        );

        ## Initiate the array of 'included components' to be passed to the next
        ## level of recursion. The first member of the array is the case where
        ## no reciprocal rearrangements are included.
        my @new_included_components = (Sanger::CGP::Rearrangement::ConnectedComponentSet->new(
            components => [
                $currently_included_components->components_array,
                $new_remaining_components->shift_first
            ]
        ));

        ## If the current component is a single rearrangement, then search and
        ## include all reciprocal rearrangements. But if the current component
        ## has multiple rearrangements (that are all linked through being
        ## balanced), then by definition the these rearrangements cannot have
        ## reciprocal real rearrangement partners other than those in the same
        ## connected components.
        my ($reciprocal_component, $cur_rg, $reciprocal_rg, $new_component_set);
        if ($remaining_components->first->size == 1) {
            ($cur_rg) = $remaining_components->first->rgs_array;
            for $reciprocal_rg ($cur_rg->reciprocal_rgs_array) {
                $reciprocal_component = $reciprocal_rg->component;
                next if $reciprocal_component->size > 1;
                die if !defined($reciprocal_component);

                ## Add the current reciprocal component only if it's not
                ## already one of the included components.
                if (
                    !$currently_included_components->contains($reciprocal_component) &&
                    !$remaining_components->contains($reciprocal_component)
                ) {
                    $new_component_set = $new_included_components[0]->copy;
                    $new_component_set->add($reciprocal_component);
                    push @new_included_components, $new_component_set;
                }
            }
        }

        return(
            map(
                {
                    make_rg_groups_by_reciprocal_rg_recursion(
                        $new_remaining_components,
                        $_,  ## RearrangementConnectedComponentSet object of included components
                        $path,
                    )
                }
                @new_included_components
            )
        );
    }
}

sub generate_rearrangement_groups {
    my $genome = shift;
    my %params = @_;
    if (!defined(blessed $genome) || blessed($genome) ne 'Sanger::CGP::CopyNumber::Segment::Genome') {
        die;
    }

    ## Algorithm: depth-first search.
    ## For each rearrangement and reciprocal rearrangement pair, find and store
    ## all paths of 'effect'. An effector path is any set of up to three
    ## rearrangements (and/or reciprocal rearrangement pairs) that are mutually
    ## located and oriented such that (one of) the rearrangement(s) of interest
    ## is potentially confounded by it.
    ##
    ## Implementation:
    ## 1. Start from low end. Find all rearrangement ends that have effect on
    ##    it.
    ## 2. Extend each path until
    ##    i: the high end is reached, or
    ##    ii: the depth of up to three is reached, at which point search is
    ##        stopped.
    ## 3. For each set of rearrangement + effector paths, collate them plus
    ##    all the combinations of reciprocal/non-reciprocal rearrangements
    ##    inclusion.
    ## 4. Make rearrangement patterns, where every balanced and/or fold-back
    ##    segment is both
    ##
    ## First find all, even redundant groups. Then remove redundant groups
    ## second.
    ##
    ## Steps 1-2 above are implemented in package 'RearrangementPath'.

    my($component, $rg, $path, $path_component, $path_component_set);
    my %rg_groups = ();
    my $new_rg_group = ();
    my $components_analysed = 0;
    for $component (@{$genome->connected_components}) {
        if ($opts{verbose} >= 2) {
            print STDERR "Number of groups is now: " . scalar(keys(%rg_groups)) . ", components analysed: " . $components_analysed++ . "\n";
            print STDERR "Rearrangement paths and groups for:\n";
            print STDERR $component->to_s;
        }

        ## First the pathless component only group.
        $new_rg_group = Sanger::CGP::Rearrangement::Group->new(
            components => [$component],
            target_component => $component,
        );
        $new_rg_group->get_normalised_rg_patterns(%opts);  ## TODO
        if (!exists($rg_groups{$new_rg_group->unique_string})) {
            $rg_groups{$new_rg_group->unique_string} = $new_rg_group;
        }

        ## Next we find all sets of components considering 'rearrangement paths'
        for $rg (values %{$component->rgs}) {
            if ($rg->is_part_of_shard_cycle(%params)) {
                next;
            }

            for $path (RearrangementPath->find_all_valid_paths(target_rg => $rg, max_depth => $opts{max_depth}, %params)) {
                ## If the current component has a single rearrangement, then we
                ## need to search for reciprocal rearrangements for the single
                ## rearrangement.
                if ($component->size == 1) {
                    $path_component_set = Sanger::CGP::Rearrangement::ConnectedComponentSet->new(components => [$component]);
                    for $path_component ($path->components_array) {
                        if ($path_component eq $component) {
                            next;
                        }
                        $path_component_set->add($path_component);
                    }
                    for $new_rg_group (make_rg_groups_by_reciprocal_rg_recursion(
                        Sanger::CGP::Rearrangement::ConnectedComponentSet->new(components => [$path_component_set->components_array]),
                        Sanger::CGP::Rearrangement::ConnectedComponentSet->new(),
                        $path,
                    )) {
                        if (!exists($rg_groups{$new_rg_group->unique_string})) {
                            $rg_groups{$new_rg_group->unique_string} = $new_rg_group;
                        }
                    }
                }
                else {
                    $path_component_set = Sanger::CGP::Rearrangement::ConnectedComponentSet->new();
                    for $path_component ($path->components_array) {
                        if ($path_component eq $component) {
                            next;
                        }
                        $path_component_set->add($path_component);
                    }
                    for $new_rg_group (make_rg_groups_by_reciprocal_rg_recursion(
                        Sanger::CGP::Rearrangement::ConnectedComponentSet->new(components => [$path_component_set->components_array]),
                        Sanger::CGP::Rearrangement::ConnectedComponentSet->new(components => [$component]),
                        $path,
                    )) {
                        if (!exists($rg_groups{$new_rg_group->unique_string})) {
                            $rg_groups{$new_rg_group->unique_string} = $new_rg_group;
                        }
                    }
                }
            }
        }
    }

    return(values %rg_groups);
}


__END__


=head1 NAME

match_rg_patterns_to_library.pl

Matching rearrangement patterns with enumerated rearrangement patterns
arising from classical rearrangement mechanisms.

=head1 SYNOPSIS

match_rg_patterns_to_library.pl [options] REARRANGEMENTS.BEDPE CN_SEGMENTS.BED

REARRANGEMENTS.BEDPE is a BEDPE file with rearrangement ID in column 7
and strands in columns 9 and 10.

CN_SEGMENTS.BEDGRAPH is a BEDGRAPH file of copy number segments.
The fourth column corresponds to the absolute copy number of the segment.
The fifth and sixth columns correspond to copy number breakpoint type
(rearrangement vs. copy number segmentation breakpoint). The seventh
column is the number of windows in the current segment.

  Basic options:
    -help                       Print this message
    -acf FLOAT                  Aberrant cell fraction [1.0]
    -ploidy FLOAT               Tumour ploidy [2.0]
    -verbose                    Print debugging messages

  Advanced options - defaults work OK
    -max_balanced_rg_dist INT   Maximum distance at which reciprocal
                                rearrangements can still be considered balanced
                                [1000]
    -max_foldback_distance INT  Maximum distance for fold-back type
                                rearrangements to be considered as purely
                                fold-back [5000].
    -min_seg_size_for_cn INT    Minimum segment size from which CNs estimates
                                are trusted and used for filtering out false
                                positive rearrangements [10000]
    -min_cn_change FLOAT        The minimum amount the copy number change
                                across a putative rearrangement call for the
                                rearrangement to be not filtered out (whenever
                                the copy number change across breakpoint can be
                                computed) [0.5]
    -keep_small_dels_and_tds    Keep TDs and deletions that are smaller than
                                min_seg_size_for_cn? [TRUE]
    -max_depth INT              Maximum length of 'rearrangement paths' [3]

