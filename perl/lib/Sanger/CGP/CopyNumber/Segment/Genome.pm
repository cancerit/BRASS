package Sanger::CGP::CopyNumber::Segment::Genome;

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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

use strict;
use warnings FATAL => 'all';
use Bio::Brass qw($VERSION);
use Sanger::CGP::CopyNumber::Segment::Array;
use Sanger::CGP::CopyNumber::Segment;
use Sanger::CGP::Rearrangement;
use Sanger::CGP::Rearrangement::ConnectedComponent;
use Sanger::CGP::Rearrangement::Group;
use List::Util qw(sum);

sub new {
    my $class = shift;
    my %params = @_;

    return bless {
        chrs => [],
        rg_of_id => ($params{rg_of_id} or undef)
    }, $class;
}

sub new_from_BEDPE_and_CN_BED {
    my $class = shift;
    my $rgs_file = shift;
    my $cn_file = shift;
    my %params = @_;
    my %rg_of_id;
    my @F;
    my(@reads_min_max_pos, @reads_clipped);

    open RGS, $rgs_file;
    while (<RGS>) {
        chomp;
        @F = split /\t/;
        @reads_min_max_pos = @reads_clipped = ();
        for (@F[12..15]) {
            /^(\d+) \((\d+)\)$/ or die;
            push @reads_min_max_pos, $1;
            push @reads_clipped, $2;
        }
        $rg_of_id{$F[6]} = Sanger::CGP::Rearrangement->new(
            low_end_dir => $F[8],
            high_end_dir => $F[9],
            id => $F[6],
            orig_data => [@F],
            reads_pos => \@reads_min_max_pos,
            reads_clip => \@reads_clipped,
        );
    }
    close RGS;

    my ($chr, $start_pos, $end_pos, $cn, $low_end_bkpt, $high_end_bkpt, $n_win);
    my (%segments_of_chr, $low_end, $high_end);
    open CN, $cn_file;
    while (<CN>) {
        chomp;
        ($chr, $start_pos, $end_pos, $cn, $low_end_bkpt, $high_end_bkpt, $n_win) = split /\t/;
        # $start_pos++;  # Convert 0-based BEDPE format start back to 1-based. # DOESN'T NEED TO BE DONE

        ## NOTE: below has to be removed later
        if ($low_end_bkpt =~ /:/ && !exists($rg_of_id{substr($low_end_bkpt, 0, -2)})) {
            warn "Rearrangement $low_end_bkpt in CN file not found in RG file. Changing into cn_bkpt";
            $low_end_bkpt = 'cn_bkpt';
        }
        if ($high_end_bkpt =~ /:/ && !exists($rg_of_id{substr($high_end_bkpt, 0, -2)})) {
            warn "Rearrangement $high_end_bkpt in CN file not found in RG file. Changing into cn_bkpt";
            $high_end_bkpt = 'cn_bkpt';
        }

        if (!(exists $segments_of_chr{$chr})) {
            $segments_of_chr{$chr} = Sanger::CGP::CopyNumber::Segment::Array->new(name => $chr);
        }

        if ($n_win eq 'NA') {
            print STDERR "Warning: encountered segment $chr:$start_pos-$end_pos (size: " . ($end_pos - $start_pos + 1) . ") with n_win 'NA'. n_win set to 1.\n" if $params{verbose} >= 2;
            $n_win = 1;
        }
        $segments_of_chr{$chr}->add_copy_number_segment(
            Sanger::CGP::CopyNumber::Segment->new(
                chr           => $segments_of_chr{$chr},
                start         => $start_pos,
                end           => $end_pos,
                cn            => ($cn eq 'NA' ? undef : $cn),
                n_win         => $n_win,
                low_end_bkpt  => $low_end_bkpt,
                high_end_bkpt => $high_end_bkpt,
                rg_of_id      => \%rg_of_id,
            )
        );
    }
    close CN;

    # Sanity check: after reading rearrangement data in, each rearrangement end
    # in %rg_of_id should be associated with a copy number segment end.
    for (keys %rg_of_id) {
        if (
            !defined($rg_of_id{$_}->low_end->{segment_end}) &&
            !defined($rg_of_id{$_}->high_end->{segment_end})
        ) {
            print STDERR "Rearrangement '$_' was not found in RG_CNS file and will be removed from \%rg_of_id.\n" if $params{verbose} >= 2;
            delete $rg_of_id{$_};
        }
        elsif (!defined($rg_of_id{$_}->low_end->{segment_end})) {
            warn "Rearrangement $_ low end was not assigned to any copy number segment!";
        }
        elsif (!defined($rg_of_id{$_}->high_end->{segment_end})) {
            warn "Rearrangement $_ high end was not assigned to any copy number segment!";
        }
    }

    return bless {
        chrs => \%segments_of_chr,
        rg_of_id => \%rg_of_id,
        connected_components => undef,
    }, $class;
}

#
# Getter methods
#
sub rg_of_id {
    my $self = shift;
    if (!defined($self->{rg_of_id})) {
        die "Attempted to call rg_of_id() to get undefined $self\->{rg_of_id}!";
    }
    return $self->{rg_of_id};
}

sub chrs {
    my $self = shift;
    # if (!defined($self->{chrs})) {
    #     die "Attempted to call $self\->chrs() to get undefined $self\->{chrs}!";
    # }
    return $self->{chrs};
}

sub chrs_array {
    my $self = shift;
    return(values %{$self->{chrs}});
}

sub connected_components {
    my $self = shift;
    if (!defined($self->{connected_components})) {
        die "Attempted to call connected_components() to get undefined $self\->{connected_components}!";
    }
    return $self->{connected_components};
}

sub connected_components_array {
    my $self = shift;
    return @{$self->connected_components};
}

sub chr_by_name {
    my $self = shift;
    my $chr_name = shift;
    if (!exists($self->chrs->{$chr_name})) {
        warn "No match for chromosome name '$chr_name' in $self\->chr_by_name()";
        return undef;
    }
    else {
        return $self->chrs->{$chr_name};
    }
}
#
# End of getter methods
#

#
# Worker subroutines
#
sub sort_segments {
    my $self = shift;
    for my $chr ($self->chrs_array) {
        $chr->sort_segments;
    }
}

sub delete_rg {
    # Deletes $rg_id from $self->{rg_of_id}, but before that makes all the
    # necessary updates in the affected copy number segment data etc.
    my $self = shift;
    my $rg_id = shift;
    my $cur_rg = $self->rg_of_id->{$rg_id};
    $cur_rg->low_end->segment_end->{bkpt} = undef;
    $cur_rg->low_end->segment_end->{boundary} = 'cn_bkpt';
    $cur_rg->low_end->segment_end_neighbour->{boundary} = 'cn_bkpt';
    $cur_rg->high_end->segment_end->{bkpt} = undef;
    $cur_rg->high_end->segment_end->{boundary} = 'cn_bkpt';
    $cur_rg->high_end->segment_end_neighbour->{boundary} = 'cn_bkpt';
    delete($self->rg_of_id->{$rg_id});
}

sub filter_rgs_by_no_cn_change {
    my $self = shift;
    my %params = @_;
    if (!exists($params{min_cn_change})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'min_cn_change'";
    }
    if (!exists($params{keep_small_dels_and_tds})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'keep_small_dels_and_tds'";
    }
    if (!exists($params{max_shard_length})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'max_shard_length'";
    }

    my($cur_rg, $has_cn_change, $both_ends_balanced, @rgs_to_delete);
    for my $rg_id (keys %{$self->rg_of_id}) {
        $cur_rg = $self->rg_of_id->{$rg_id};

        if ($params{keep_small_dels_and_tds} && ($cur_rg->is_small_td(%params) || $cur_rg->is_small_del(%params))) {
            next;
        }

        # If the segments are too short for copy number to be estimated properly
        if (
            !$cur_rg->is_foldback(%params) &&
            (
                (
                    $cur_rg->low_end->segment->length < 1000 ||
                    $cur_rg->low_end->segment_end_neighbour->segment->length < 1000
                ) &&
                (
                    $cur_rg->high_end->segment->length < 1000 ||
                    $cur_rg->high_end->segment_end_neighbour->segment->length < 1000
                )
            )
        ) {
            $has_cn_change = 1;
        }
        elsif (defined($cur_rg->weighted_avg_cn_change_across_rg(%params))) {
            $has_cn_change = $cur_rg->weighted_avg_cn_change_across_rg(%params) >= $params{min_cn_change};
        }
        else {
            ## The fact that we got here means that the current rearrangement
            ## is a fold-back but had the exact same breakpoint with some other
            ## rearrangments.
            $has_cn_change = 1;
        }

        if (!$cur_rg->is_foldback(%params)) {
            $both_ends_balanced = defined($cur_rg->low_end->balanced_bkpt_partner_rg_end(%params)) &&
                                  defined($cur_rg->high_end->balanced_bkpt_partner_rg_end(%params));

            if (!$has_cn_change && !$both_ends_balanced) {
                push @rgs_to_delete, $rg_id;
            }
        }
        elsif (
            !$has_cn_change &&
            !($cur_rg->low_end->dir eq '+' && defined($cur_rg->high_end->balanced_bkpt_partner_rg_end(%params))) &&
            !($cur_rg->low_end->dir eq '-' && defined($cur_rg->low_end->balanced_bkpt_partner_rg_end(%params)))
        ) {
            push @rgs_to_delete, $rg_id;
        }
    }

    for my $rg_id (@rgs_to_delete) {
        $self->delete_rg($rg_id);
    }

    return scalar(@rgs_to_delete);
}

sub remove_shard_sequence_bypassing_rgs {
    my $self = shift;
    my %params = @_;
    my @rgs_to_delete = ();
    my $rg_id;
    for $rg_id (keys %{$self->rg_of_id}) {
        if ($self->rg_of_id->{$rg_id}->is_shard_bypassing(%params)) {
            push @rgs_to_delete, $rg_id;
        }
    }

    for $rg_id (@rgs_to_delete) {
        $self->delete_rg($rg_id);
    }

    return scalar(@rgs_to_delete);
}

sub remove_small_cn_bkpt_segments {
    my $self = shift;
    my %params = @_;
    for my $chr ($self->chrs_array) {
        $chr->remove_small_cn_bkpt_segments(%params);
    }
}

sub remove_cn_bkpts_without_cn_change {
    my $self = shift;
    my %params = @_;
    for my $chr ($self->chrs_array) {
        $chr->remove_cn_bkpts_without_cn_change(%params);
    }
}

sub preprocess_rearrangement_and_copy_number_data {
    my $self = shift;
    my %params = @_;
    $self->sort_segments;

    my $changed = 1;
    my $iterations = 0;
    while ($changed) {
        $changed = 0;
        $iterations++;
        if ($params{verbose} >= 2) {
            print STDERR "  Iterations: $iterations...\n";
        }
        $changed += $self->remove_cn_bkpts_without_cn_change(%params);
        $changed += $self->remove_small_cn_bkpt_segments(%params);
        # print STDERR "Removing shard bypassing rearrangements has been disabled...\n";
        $changed += $self->remove_shard_sequence_bypassing_rgs(%params);
        $changed += $self->filter_rgs_by_no_cn_change(%params);
    }

    return $iterations;
}

# Function for marking balanced rearrangement pairs.
sub mark_reciprocal_rearrangements {
    my $self = shift;
    my %params = @_;
    my @rearrangements = values %{$self->rg_of_id};
    my $i;
    my @reciprocal_rg_ends;
    for $i (0..($#rearrangements-1)) {
        ## Start the reciprocal rearrangement search from the low end.
        if ($rearrangements[$i]->low_end->is_fwd) {
            @reciprocal_rg_ends = $rearrangements[$i]->low_end->get_us_effector_rg_ends(%params);
        }
        else {
            @reciprocal_rg_ends = $rearrangements[$i]->low_end->get_ds_effector_rg_ends(%params);
        }
        for (@reciprocal_rg_ends) {
            if ($rearrangements[$i]->id ge $_->id) {
                next;
            }
            elsif ($rearrangements[$i]->is_reciprocal_with($_->rg)) {
                $rearrangements[$i]->add_reciprocal_rgs($_->rg);
                $_->rg->add_reciprocal_rgs($rearrangements[$i]);
            }
        }
    }
}

sub compute_connected_components {
    ## A breadth-first search algorithm for connected components of
    ## rearrangements. Currently rearrangements become connected only
    ## through being balanced.
    my $self = shift;
    my %params = @_;

    my %unseen_rearrangements = %{$self->rg_of_id};
    my $current_component;
    my @connected_components = ();  ## A list of all connected components.
    my @find_neighbours_of = ();    ## A helper list for keeping tab of unanalysed rearrangements
                                    ## in the current component. FIFO.
    my ($cur_rg_id, $cur_rg, $partner_rg_end);

    for my $rg_id (keys %unseen_rearrangements) {
        # Has this been already 'seen' in previous iterations?
        if (!exists($unseen_rearrangements{$rg_id})) {
            next;
        }
        $cur_rg = $unseen_rearrangements{$rg_id};

        # Start breadth-first search with $rg_id. First initiate a new 'current component'.
        $current_component = Sanger::CGP::Rearrangement::ConnectedComponent->new(
            rgs => { $rg_id => $cur_rg }
        );
        die if defined($cur_rg->component);
        $cur_rg->set_component($current_component);
        @find_neighbours_of = ($rg_id);
        delete $unseen_rearrangements{$rg_id};

        # Then iterate through the breadth-first search.
        while (@find_neighbours_of) {
            $cur_rg_id = shift @find_neighbours_of;
            $cur_rg    = $current_component->rgs->{$cur_rg_id};

            # Find potential neighbours for low end...
            $partner_rg_end = $cur_rg->low_end->balanced_bkpt_partner_rg_end(%params);
            if (
                defined($partner_rg_end) &&
                exists($unseen_rearrangements{$partner_rg_end->id})
            ) {
                push @find_neighbours_of, $partner_rg_end->id;
                delete $unseen_rearrangements{$partner_rg_end->id};
                $current_component->add_rgs($partner_rg_end->rg);
                die if defined($partner_rg_end->rg->component);
                $partner_rg_end->rg->set_component($current_component);
            }

            # ... and high end.
            $partner_rg_end = $cur_rg->high_end->balanced_bkpt_partner_rg_end(%params);
            if (
                defined($partner_rg_end) &&
                exists($unseen_rearrangements{$partner_rg_end->id})
            ) {
                push @find_neighbours_of, $partner_rg_end->id;
                delete $unseen_rearrangements{$partner_rg_end->id};
                $current_component->add_rgs($partner_rg_end->rg);
                die if defined($partner_rg_end->rg->component);
                $partner_rg_end->rg->set_component($current_component);
            }
        }

        # Store the current component and prepare for the next available rearrangement.
        push @connected_components, $current_component;
    }

    $self->{connected_components} = \@connected_components;
}

sub print_rearrangements_in_bedpe {
    my ($self, $fh, %params) = @_;
    $fh = *STDOUT unless(defined $fh);
    my $rg;
    for (keys %{$self->rg_of_id}) {
        $rg = $self->rg_of_id->{$_};
        print $fh join(
            "\t",
            @{$rg->orig_data}[0..5],
            $rg->id,
            $rg->orig_data->[7],
            $rg->low_end->dir,
            $rg->high_end->dir,
            ($rg->low_end->cn_across_bkpt(%params) || 'NA'),
            ($rg->high_end->cn_across_bkpt(%params) || 'NA'),
            @{$rg->orig_data}[10..$#{$rg->orig_data}],
        ) . "\n";
    }
}

sub print_rg_cns_bedpe {
    my ($self, $fh) = @_;
    for my $chr (sort {$a->name cmp $b->name} values($self->chrs)) {
        $chr->print_rg_cns_bedpe($fh);
    }
}

sub print_naive_classifications {
    my $self = shift;
    my %params = @_;
    my ($c, $classification, $event_id);
    $event_id = 0;

    my %components_in_complex_regions = ();
    my($event, $component);
    for $event (@{$self->{complex_events}}) {  ##TODO
        for $component ($event->components_array) {
            for ($component->rgs_array) {
                print join(
                    "\t",
                    $_->low_end->chr_name,
                    $_->low_end->pos - 1,
                    $_->low_end->pos,
                    $_->high_end->chr_name,
                    $_->high_end->pos - 1,
                    $_->high_end->pos,
                    ("$event_id:" . $_->id . ':chromothripsis'),
                    1,
                    $_->low_end->dir,
                    $_->high_end->dir,
                ) . "\n";
            }

            $components_in_complex_regions{$component} = 1;
        }
        $event_id++;
    }

    for my $c ($self->connected_components_array) {
        if (exists($components_in_complex_regions{$c})) {
            next;
        }

        $classification = $c->naive_classification(%params);
        for ($c->rgs_array) {
            print join(
                "\t",
                $_->low_end->chr_name,
                $_->low_end->pos - 1,
                $_->low_end->pos,
                $_->high_end->chr_name,
                $_->high_end->pos - 1,
                $_->high_end->pos,
                ("$event_id:" . $_->id . ":$classification"),
                1,
                $_->low_end->dir,
                $_->high_end->dir
            ) . "\n";
        }
        $event_id++;
    }
}

sub print_unique_breakpoints {
    my $self = shift;
    my %params = @_;
    for my $chr (sort { $a->name cmp $b->name } $self->chrs_array) {
        $chr->print_unique_breakpoints(%params);
    }
}

sub find_rg_pattern_motifs {
    my $self = shift;
    my %params = @_;
    if (!exists($params{within})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'within'";
    }
    if (!exists($params{away_from})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'away_from'";
    }
    if (!exists($params{size})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'size'";
    }

    # Strategy:
    # Store found motifs in a hash.
    # Then go through each rearrangement in turn.
    # Find neighbours of each rearrangment, and grow the size of each group
    # until the desired size is found.
    my %selected_groups = ();
    my %selected_rgs;
    my $chr;
    my $seg;
    my @rg_ends;
    my $group;
    my($rg_end, $rg);
    my $have_non_involved_rgs_nearby;
    my @neighbour_candidates;

    sub rg_neighbours {
        my $rg = shift;
        my %params = shift;
        my %neighbours;
        for ($rg->low_end->closest_rg_end(%params), $rg->high_end->closest_rg_end(%params)) {
            $neighbours{$_->rg->id} = $_->rg;
        }
        return values(%neighbours);
    }

    sub iterate_further_rgs {
        my %included_rgs = @{shift()};
        my $next_rg_to_be_included = shift;
        my $depth = shift;
        my %params = @_;
        my $rg;

        $included_rgs{$next_rg_to_be_included->id} = $next_rg_to_be_included;
        if ($depth == 1) {
            # Check if any rearrangement end of the current included rearrangements are
            # too close to a non-included rearrangement.
            for $rg (values %included_rgs) {
                if (
                    defined($rg->low_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs])) ||
                    defined($rg->high_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs]))
                ) {
                    return();
                }
            }

            return(
                join(',', sort(keys %included_rgs)),
                Sanger::CGP::Rearrangement::Group->new(
                    components => [map { $_->component } values(%included_rgs)]
                )
            );
        }
        else {
            # Find new potential neighbours and iterate further.
            my @iterated_groups = ();
            my %neighbours = ();

            for $rg (values %included_rgs) {
                for (
                    $rg->low_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs]),
                    $rg->high_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs])
                ) {
                    if (defined($_)) {
                        $neighbours{$_->rg->id} = $_->rg;
                    }
                }
            }

            for $rg (values %neighbours) {
                push @iterated_groups, iterate_further_rgs([%included_rgs], $rg, $depth-1, %params);
            }

            return @iterated_groups;
        }
    }

    for (values $self->rg_of_id) {
        %selected_groups = (%selected_groups, iterate_further_rgs([], $_, $params{size}, %params));
    }

#     for $chr ($self->chrs_array) {
#         @rg_ends = ();
#         for $seg ($chr->segments_array) {
#             if (defined($seg->low_end->bkpt)) {
#                 push @rg_ends, $seg->low_end->bkpt;
#             }
#             if (defined($seg->high_end->bkpt)) {
#                 push @rg_ends, $seg->high_end->bkpt;
#             }
#         }
#
#         # Strategy for each encountered rearrangement end:
#         # 1. Collect an array of $params{size} successive rearrangements
#         # 2. Check that none of the involved rearrangements ends are within
#         #    $params{away_from} base pairs of a non-involved rearrangement end.
#         # 3. Make sure rearrangements are within a window of $params{within} size
#         my $next_idx;
#         GROUP_CANDIDATE:
#         for (0..($#rg_ends - $params{size})) {
#             %selected_rgs = ();
#             $next_idx = $_-1;
#             while (scalar(keys %selected_rgs) < $params{size} && $next_idx+1 < $#rg_ends) {
#                 $next_idx++;
#                 if ($rg_ends[$next_idx]->segment_end->pos - $rg_ends[$_]->segment_end->pos > $params{within}) {
#                     next GROUP_CANDIDATE;
#                 }
#                 $selected_rgs{$rg_ends[$next_idx]->rg->id} = $rg_ends[$next_idx]->rg;
#             }
#             if (scalar(keys %selected_rgs) < $params{size}) {
#                 next;
#             }
#
#             # Check if the current group is already handled/found
#             for $rg_end (@rg_ends[$_..($_+$params{size}-1)]) {
#                 $selected_rgs{$rg_end->rg->id} = $rg_end->rg;
#             }
#             if (exists($selected_groups{ join(',', sort(keys %selected_rgs)) })) {
#                 next;
#             }
#
#             # Make sure none of the involved rearrangements are within near
#             # non-involved rearrangement ends.
#             $have_non_involved_rgs_nearby = 0;
#             for $rg (values %selected_rgs) {
#                 if (
#                     defined($rg->low_end->closest_rg_end(within => $params{away_from}, exclude_rgs => [keys %selected_rgs])) ||
#                     defined($rg->high_end->closest_rg_end(within => $params{away_from}, exclude_rgs => [keys %selected_rgs]))
#                 ) {
#                     $have_non_involved_rgs_nearby = 1;
#                     last;
#                 }
#             }
#             if ($have_non_involved_rgs_nearby) {
#                 next;
#             }
#
#             # Finally store the group
#             $selected_groups{join(',', sort(keys %selected_rgs))} = RearrangementGroup->new(
#                 components => [map { $_->component } values(%selected_rgs)]
#             );
#
#             # Sanity check - should never happen if $params{within} < $params{away_from}
#             if (scalar($selected_groups{join(',', sort(keys %selected_rgs))}->rgs_array) > $params{size}) {
#                 die;
#             }
#         }
#     }
#

    map { $_->get_normalised_rg_patterns(%params) } values(%selected_groups);
    return values(%selected_groups);
}
#
# End of worker subroutines
#


1;
