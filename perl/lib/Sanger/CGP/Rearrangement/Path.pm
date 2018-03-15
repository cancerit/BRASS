package Sanger::CGP::Rearrangement::Path;

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

use warnings FATAL => 'all';
use strict;
use Bio::Brass qw($VERSION);

use Scalar::Util qw(blessed);

sub find_all_valid_paths {
    my $class = shift;
    my %params = @_;
    if (!exists($params{target_rg}) || !exists($params{max_depth})) {
        die;
    }

    # Generate the first RearrangementPath object
    my $path = bless {
        target_rg => $params{target_rg},
        path => [],
    }, $class;

    my $target_rg_ends = {};
    for ($params{target_rg}->high_end->get_us_effector_rg_ends(%params), $params{target_rg}->high_end->get_ds_effector_rg_ends(%params)) {
        $target_rg_ends->{$_} = 1;
    }

    ## The path discovery is carried out below.
    return find_paths(
        path => $path,
        in_cis => 0,
        depth => 0,
        max_depth => $params{max_depth},
        target_rg_ends => $target_rg_ends,
        min_seg_size_for_cn => $params{min_seg_size_for_cn},
        max_shard_length => $params{max_shard_length},
    );
}

#
# Getter and setter methods
#
sub target_rg {
    my $self = shift;
    return $self->{target_rg};
}

sub path {
    my $self = shift;
    return $self->{path};
}

sub components_array {
    ## The components which include all the rearrangements as part of the path.
    my $self = shift;
    my %components = ();
    for my $rg_end (@{$self->path}) {
        if (!exists($components{$rg_end->rg->component})) {
            $components{$rg_end->rg->component} = $rg_end->rg->component;
        }
    }
    return values(%components);
}

sub length {
    my $self = shift;
    return scalar(@{$self->path});
}

## Below is not used??
# sub extend_path {
#     my $self = shift;
#     my $target_rg_end = shift;
#     if (!defined(blessed $target_rg_end) || blessed($target_rg_end) ne "RearrangementEnd") {
#         die "Argument of type 'RearrangementEnd' is needed for $self\->extend_path()";
#     }
#
#     push @{$self->{path}}, $target_rg_end;
# }

sub last_node {
    my $self = shift;
    if (scalar(@{$self->path}) < 1) {
        return undef;
    }
    else {
        return $self->path->[$#{$self->path}];
    }
}

sub to_s {
    my $self = shift;
    my %params = @_;
    my $output = 'Target: ' . $self->target_rg->to_s . "\n";
    my ($mate, $shard_count);
    for (@{$self->path}) {
        ($mate, $shard_count) = $_->traverse_across_shards(%params);
        $output .= $_->to_s . "\t" . $_->id . "\t$shard_count shards\n\t" . $mate->to_s . "\n";
    }

    return $output;
}

sub print {
    my $self = shift;
    my %params = @_;
    print $self->to_s(%params);
}
sub p { $_[0]->print; }
#
# End of getter and setter methods
#

#
# Helper methods
#
sub extend_and_clone {
    my $self = shift;
    my $next_rg_end = shift;
    if (!defined(Scalar::Util::blessed $next_rg_end) || Scalar::Util::blessed($next_rg_end) ne 'Sanger::CGP::Rearrangement::End') {
        die "Argument of type 'Sanger::CGP::Rearrangement::End' is needed for $self\->extend_and_clone()";
    }

    my $new_path = bless {
        target_rg => $self->target_rg,
        path => [@{$self->path}, $next_rg_end],
    }, ref($self);
    return $new_path;
}

sub check_is_path_complete {
    my $self = shift;
    my %params = @_;
    if (!@{$self->path}) {
        return 0;
    }
    my $source_rg_end = $self->target_rg->low_end;
    my $target_rg_end = $self->target_rg->high_end;

    ## Check whether the first node of the path is oriented properly with
    ## respect to the low end of the target rearrangement.
    if ($self->path->[0]->chr ne $source_rg_end->chr) {
        die;
    }
    if (
        !($self->path->[0]->pos <= $source_rg_end->pos && $self->path->[0]->is_rev) &&
        !($self->path->[0]->pos >= $source_rg_end->pos && $self->path->[0]->is_fwd)
    ) {
        die;
    }

    my $mate;
    for (0..($#{$self->path}-1)) {
        ## Check whether the cis-requirement of each transition in
        ## the path is respected.
        ($mate) = $self->path->[$_]->traverse_across_shards(%params);
        if (!$mate->can_be_in_cis_with($self->path->[$_+1])) {
            die;
        }
    }

    ## Check if the final node of the path is oriented properly with
    ## respect to the high end of the target rearrangement.
    my $final_rg_end = ($self->last_node->traverse_across_shards(%params))[0];
    if ($final_rg_end->chr ne $target_rg_end->chr) {
        die;
    }
    if (
        !($final_rg_end->pos <= $target_rg_end->pos && $final_rg_end->is_rev) &&
        !($final_rg_end->pos >= $target_rg_end->pos && $final_rg_end->is_fwd)
    ) {
        die;
    }
}

sub alters_rg_end_relative_positioning {
    my $self = shift;
    my %params = @_;
    my $path = $self->path;
    my $low_end = $self->target_rg->low_end;
    my $high_end = $self->target_rg->high_end;

    ## Otherwise does the path potentially lead to another chromosome?
    if ($low_end->chr_name ne $high_end->chr_name) {
        return 1;
    }

    ## Check whether order or orientation of the two rearrangements is altered.
    if ($high_end->pos >= $low_end->pos) {
        ## High end seemingly downstream of low end but actually could be upstream.
        if ($path->[0]->pos < $low_end->pos) {
            return 1;
        }
    }
    else {
        ## High end seemingly upstream of low end but actually could be downstream.
        ## THIS SHOULD NEVER HAPPEN!!!
        # if ($path->[0]->pos > $low_end->pos) {
        #     return 1;
        # }
        die;
    }

    ## Otherwise is there a net change in orientation?
    if ($self->path->[0]->dir eq ($self->last_node->traverse_across_shards(%params))[0]->dir) {
        return 1;
    }


    return 0;
}

## Below function recursively extends the paths in order to find alternative
## 'paths' from the low end of the target rearrangement to the high end.
sub find_paths {
    my %params = @_;
    my $in_cis = $params{in_cis};  ## Whether the next step in the path has to be in cis
                                   ## with the previous step of the path.
    my $cur_path = $params{path};
    my $depth = $params{depth} + 1;
    my $max_depth = $params{max_depth};
    my $target_rg_ends = $params{target_rg_ends};
    my @effector_rg_ends;
    my @output_paths = ();

    ## The current 'anchor' RearrangementEnd against which we try to find
    ## effector rearrangement ends.
    my $target_rg_end;
    if (!defined($cur_path->last_node)) {
        $target_rg_end = $cur_path->target_rg->low_end;
    }
    else {
        ($target_rg_end) = $cur_path->last_node->traverse_across_shards(%params);
    }

    if (!$in_cis) {
        @effector_rg_ends = (
            $target_rg_end->get_us_effector_rg_ends(%params),
            $target_rg_end->get_ds_effector_rg_ends(%params)
        );
    }
    elsif ($target_rg_end->is_fwd) {
        @effector_rg_ends = $target_rg_end->get_us_effector_rg_ends(%params);
    }
    else {
        @effector_rg_ends = $target_rg_end->get_ds_effector_rg_ends(%params);
    }

    ## For each extended path, see if we were able to reach any of the targets.
    my $new_path;
    for my $eff (@effector_rg_ends) {
        $new_path = $cur_path->extend_and_clone($eff);
        if (exists($target_rg_ends->{($eff->traverse_across_shards(%params))[0]})) {
            $new_path->check_is_path_complete(%params);  ## A sanity check
            if ($new_path->alters_rg_end_relative_positioning(%params)) {
                push @output_paths, $new_path;
            }
        }

        if ($depth < $max_depth) {
            push @output_paths, find_paths(
                path => $new_path,
                in_cis => 1,
                depth => $depth,
                max_depth => $max_depth,
                target_rg_ends => $target_rg_ends,
                min_seg_size_for_cn => $params{min_seg_size_for_cn},
                max_shard_length => $params{max_shard_length},
            );
        }
    }

    return @output_paths;
}
#
# End of helper methods
#

1;
