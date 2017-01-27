package Sanger::CGP::Rearrangement::End;

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
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

use warnings FATAL => 'all';
use strict;
use Scalar::Util qw(refaddr blessed);
use List::Util qw(min max);
use Bio::Brass qw($VERSION);

sub new {
    my $class = shift;
    my %params = @_;
    return bless {
        end            => $params{end},  # 'low' or 'high'
        dir            => $params{dir},
        segment_end    => ($params{segment_end} || undef),
        rg             => ($params{rg} || undef),
        min_reads_pos  => ($params{min_reads_pos} || undef),
        min_reads_clip => ($params{min_reads_clipped} || undef),
        max_reads_pos  => ($params{max_reads_pos} || undef),
        max_reads_clip => ($params{max_reads_clipped} || undef),
    }, $class;
}

#
# Getter methods
#
sub end {
    my $self = shift;
    if (!defined($self->{end})) {
        die "Attempted to call $self\->end() to get undefined $self\->{end}!";
    }
    if ($self->{end} ne 'low' && $self->{end} ne 'high') {
        die "Attempted to call $self\->end() with the value $self\->{end} eq '$self->{end}' - expected 'low' or 'high'";
    }
    return $self->{end};
}

sub dir {
    my $self = shift;
    if (!defined($self->{dir})) {
        die "Attempted to call dir() to get undefined $self\->{dir}!";
    }
    if ($self->{dir} ne '+' && $self->{dir} ne '-') {
        die "Attempted to call $self\->dir() but $self\->{dir} is neither '+' nor '-'!";
    }
    return $self->{dir};
}

sub is_l {
    my $self = shift;
    return $self->dir eq 'low';
}

sub is_h {
    my $self = shift;
    return $self->dir eq 'high';
}

sub is_fwd {
    my $self = shift;
    if ($self->dir ne '+' && $self->dir ne '-') {
        die "Called $self\->is_fwd() with $self\->{dir} not in ('+', '-')";
    }
    return $self->dir eq '+';
}

sub is_rev {
    my $self = shift;
    if ($self->dir ne '+' && $self->dir ne '-') {
        die "Called $self\->is_fwd() with $self\->{dir} not in ('+', '-')";
    }
    return $self->dir eq '-';
}

sub segment_end {
    my $self = shift;
    if (!defined($self->{segment_end})) {
        die "Attempted to call segment_end() to get undefined $self\->{segment_end}!";
    }
    return $self->{segment_end};
}

sub rg {
    my $self = shift;
    if (!defined($self->{rg})) {
        die "Attempted to call rg() to get undefined $self\->{rg}!";
    }
    return $self->{rg};
}

sub min_reads_pos {
    my $self = shift;
    if (!defined($self->{min_reads_pos})) {
        die "Attempted to call min_reads_pos() to get undefined $self\->{min_reads_pos}!";
    }
    return $self->{min_reads_pos};
}

sub min_reads_clip {
    my $self = shift;
    if (!defined($self->{min_reads_clip})) {
        die "Attempted to call min_reads_clip() to get undefined $self\->{min_reads_clip}!";
    }
    return $self->{min_reads_clip};
}

sub max_reads_pos {
    my $self = shift;
    if (!defined($self->{max_reads_pos})) {
        die "Attempted to call max_reads_pos() to get undefined $self\->{max_reads_pos}!";
    }
    return $self->{max_reads_pos};
}

sub max_reads_clip {
    my $self = shift;
    if (!defined($self->{max_reads_clip})) {
        die "Attempted to call max_reads_clip() to get undefined $self\->{max_reads_clip}!";
    }
    return $self->{max_reads_clip};
}
#
# End of getter methods
#


#
# Helper subroutines
#
sub eq {
    my $self = shift;
    my $target = shift;
    return Scalar::Util::refaddr($self) == Scalar::Util::refaddr($target);
}

sub segment {
    my $self = shift;
    return $self->segment_end->segment;
}

sub chr {
    my $self = shift;
    return $self->segment->chr;
}

sub chr_name {
    my $self = shift;
    return $self->chr->name;
}

sub pos {
    my $self = shift;
    return $self->segment_end->pos;
}

sub segment_end_neighbour {
    my $self = shift;
    if ($self->segment_end->end eq 'high') {
        return $self->segment->next_seg->low_end;
    }
    elsif ($self->segment_end->end eq 'low') {
        return $self->segment->prev_seg->high_end;
    }
    else {
        die;
    }
}

sub id {
    my $self = shift;
    return $self->rg->id;
}

sub mate {
    my $self = shift;
    if ($self->end eq 'low') {
        return $self->rg->high_end;
    }
    else {
        return $self->rg->low_end;
    }
}

sub to_s {
    my $self = shift;
    return $self->chr_name . ':' . $self->pos . ':' . $self->dir;
}

sub print {
    my $self = shift;
    print $self->to_s . "\n";
}

sub is_foldback {
    my $self = shift;
    my %params = @_;
    return $self->rg->is_foldback(%params);
}

sub is_reciprocal_with {
    my $self = shift;
    my $target = shift;
    if (!defined($target) || Scalar::Util::blessed($target) ne 'Sanger::CGP::Rearrangement::End') {
        die "A argument of type 'Sanger::CGP::Rearrangement::End' is needed for $self\->is_reciprocal_with()";
    }

    if ($self->chr ne $target->chr) {
        return 0;
    }

    if ($self->is_fwd && $target->is_rev && $self->pos < $target->pos) {
        return 1;
    }

    if ($self->is_rev && $target->is_fwd && $target->pos < $self->pos) {
        return 1;
    }

    return 0;
}

sub is_on_shard {
    my $self = shift;
    my %params = @_;
    return $self->segment->is_shard(%params);
}

sub shard_partner {
    my $self = shift;
    my %params = @_;
    if (!$self->is_on_shard(%params)) {
        die "Attempted to call $self\->shard_partner() on a $self that is not on a shard";
    }

    if ($self->dir eq '+') {
        return $self->segment->low_end->bkpt;
    }
    else {
        return $self->segment->high_end->bkpt;
    }
}

sub traverse_across_shards {
    my $self = shift;
    my %params = @_;
    my $cur_rg_end = $self;
    my $shard_count = 0;
    while ($cur_rg_end->mate->is_on_shard(%params)) {
        $cur_rg_end = $cur_rg_end->mate->shard_partner(%params);
        $shard_count++;

        if ($cur_rg_end == $self) {
            ## If we get here then we got a cyclical shard sequence
            return(undef, $shard_count);
        }
    }

    return($cur_rg_end->mate, $shard_count);
}

sub can_be_in_cis_with {
    my $self = shift;
    my $target = shift;
    if (!defined($target) || Scalar::Util::blessed($target) ne 'Sanger::CGP::Rearrangement::End') {
        die "A argument of type 'Sanger::CGP::Rearrangement::End' is needed for $self\->can_be_in_cis_with()";
    }

    return(
        $self->chr eq $target->chr &&
        (
            ( $self->pos <= $target->pos && $self->is_rev && $target->is_fwd ) ||
            ( $self->pos >= $target->pos && $self->is_fwd && $target->is_rev )
        )
    );
}
#
# End of helper subroutines
#

#
# Worker subroutines
#
sub rg_side_cn {
    my $self = shift;
    my %params = @_;
    if (!$params{min_seg_size_for_cn}) {
        die "Parameter 'min_seg_size_for_cn' is needed for $self\-rg_side_cn()";
    }

    # if (
    #     !$self->rg->is_small_td(%params) &&
    #     $self->segment->length < $params{min_seg_size_for_cn}
    # ) {
    #     return undef;
    # }
    # return $self->segment->cn;

    if ($self->segment->is_bal_rg_overlap(%params)) {
        if ($self->is_fwd) {
            return $self->segment->prev_seg->cn;
        }
        else {
            return $self->segment->next_seg->cn;
        }
    }
    else {
        return $self->segment->cn;
    }
}

sub non_rg_side_cn {
    my $self = shift;
    my %params = @_;
    if (!$params{min_seg_size_for_cn}) {
        die "Parameter 'min_seg_size_for_cn' is needed for $self\-non_rg_side_cn()";
    }

    # if (
    #     !$self->rg->is_small_del(%params) &&
    #     $self->segment_end_neighbour->segment->length < $params{min_seg_size_for_cn}
    # ) {
    #     return undef;
    # }
    # return $self->segment_end_neighbour->segment->cn;

#     ## Handle balanced rearrangement ends slightly differently
#     if ($self->has_balanced_bkpt_partner_rg_end) {
#         if ($self->segment_end_neighbour->next_seg->segment->length < $params{min_seg_size_for_cn}) {
#             return undef;
#         }
#         return $self->segment_end_neighbour->next_seg->segment->cn;
#     }
#     else {
#         if ($self->segment_end_neighbour->segment->length < $params{min_seg_size_for_cn}) {
#             return undef;
#         }
#         return $self->segment_end_neighbour->segment->cn;
#     }

    # Doesn't matter if the current rg end is an overlap balanced rearrangement or not.
    if ($self->is_fwd) {
        return $self->segment->next_seg->cn;
    }
    else {
        return $self->segment->prev_seg->cn;
    }
}

sub cn_across_bkpt {
    my $self = shift;
    my %params = @_;
    my $rg_side_cn = $self->rg_side_cn(%params);
    my $non_rg_side_cn = $self->non_rg_side_cn(%params);
    if (!defined($rg_side_cn) || !defined($non_rg_side_cn)) {
        return undef;
    }

    return $rg_side_cn - $non_rg_side_cn;
}

sub relative_rg_side_cn_var {
    my $self = shift;
    my %params = @_;
    if (!$params{ploidy}) {
        die "Parameter 'ploidy' is needed for $self\-cn_change_across_rg()";
    }
    if (!$params{acf}) {
        die "Parameter 'acf' is needed for $self\-cn_change_across_rg()";
    }

    my $rg_side_cn = $self->rg_side_cn(%params);
    if (!defined($rg_side_cn)) {
        die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
        return undef;
    }
    $rg_side_cn = 0 if $rg_side_cn < 0;
    my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                              ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );

    return $adjusted_rg_side_cn / $self->segment->n_win;
}

sub relative_cn_diff_var {
    my $self = shift;
    my %params = @_;
    if (!$params{ploidy}) {
        die "Parameter 'ploidy' is needed for $self\-cn_change_across_rg()";
    }
    if (!$params{acf}) {
        die "Parameter 'acf' is needed for $self\-cn_change_across_rg()";
    }

    if ($self->segment->is_bal_rg_overlap(%params)) {
        # We can't estimate the copy number because the breakpoint
        # is completely balanced (if not even more).
        return 1e9;
    }

    ## Fold-back rearrangements must be treated a bit differently
    if ($self->rg->is_foldback(%params)) {
        my $rg_side_cn;
        my $non_rg_side_cn;

        if ($self->is_fwd) {
            $rg_side_cn = $self->rg->low_end->rg_side_cn(%params);
        }
        else {
            $rg_side_cn = $self->rg->high_end->rg_side_cn(%params);
        }
        if (!defined($rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
            return undef;
        }
        $rg_side_cn = 0 if $rg_side_cn < 0;
        if ($self->is_fwd) {
            $non_rg_side_cn = $self->rg->high_end->non_rg_side_cn(%params);
        }
        else {
            $non_rg_side_cn = $self->rg->low_end->non_rg_side_cn(%params);
        }
        if (!defined($non_rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->non_rg_side_cn";
            return undef;
        }
        $non_rg_side_cn = 0 if $non_rg_side_cn < 0;

        my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                                  ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );
        my $adjusted_non_rg_side_cn = ( 2*(1-$params{acf}) + $non_rg_side_cn*$params{acf} ) /
                                      ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );

        # Adjusted the expected read count variances by the length of the segments.
        my $adjusted_cn_diff_var = 0;
        if ($self->is_fwd) {
            $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->rg->low_end->segment->n_win +
                                    $adjusted_non_rg_side_cn / $self->rg->high_end->segment_end_neighbour->segment->n_win;
        }
        else {
            $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->rg->high_end->segment->n_win +
                                    $adjusted_non_rg_side_cn / $self->rg->low_end->segment_end_neighbour->segment->n_win;
        }

        return $adjusted_cn_diff_var;
    }
    else {
        ## The adjusted copy numbers are proportional to the expected read counts,
        ## which is proportional to the variance of the expected read counts.
        my $rg_side_cn = $self->rg_side_cn(%params);
        if (!defined($rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
            return undef;
        }
        $rg_side_cn = 0 if $rg_side_cn < 0;
        my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                                  ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );
        my $non_rg_side_cn = $self->non_rg_side_cn(%params);
        if (!defined($non_rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->non_rg_side_cn";
            return undef;
        }
        $non_rg_side_cn = 0 if $non_rg_side_cn < 0;
        my $adjusted_non_rg_side_cn = ( 2*(1-$params{acf}) + $non_rg_side_cn*$params{acf} ) /
                                      ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );

        # Adjusted the expected read count variances by the length of the segments.
        my $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->segment->n_win +
                                   $adjusted_non_rg_side_cn / $self->segment_end_neighbour->segment->n_win;
        return $adjusted_cn_diff_var;
    }
}

sub relative_cn_diff_stdev {
    my $self = shift;
    my %params = @_;
    return sqrt($self->relative_cn_diff_var(%params));
}

sub balanced_bkpt_partner_rg_end {
    # Two rearrangements are defined as balanced if they
    # are at the same copy number breakpoint pointing
    # towards each other.
    my $self = shift;
    my %params = @_;
    if ($self->segment->is_bal_rg_overlap(%params)) {
        if ($self->is_fwd) {
            return $self->segment->low_end->bkpt;
        }
        else {
            return $self->segment->high_end->bkpt;
        }
    }
    elsif ($self->is_fwd && $self->segment->next_seg->is_bal_rg_gap(%params)) {
        return $self->segment->next_seg->next_seg->low_end->bkpt;
    }
    elsif ($self->is_rev && $self->segment->prev_seg->is_bal_rg_gap(%params)) {
        return $self->segment->prev_seg->prev_seg->high_end->bkpt;
    }
    else {
        return undef;
    }
}

sub has_balanced_bkpt_partner_rg_end {
    my $self = shift;
    my %params = @_;
    return defined($self->balanced_bkpt_partner_rg_end(%params));
}

sub get_us_effector_rg_ends {
    my $self = shift;
    my %params = @_;
    my @effector_rg_ends = ();

    if (
        $self->is_fwd &&
        defined($self->segment->low_end->bkpt) &&
        !$self->mate->eq($self->segment->low_end->bkpt)
    ) {
        push @effector_rg_ends, $self->segment->low_end->bkpt;
    }

    my $cur_seg = $self->segment;
    while (defined($cur_seg->prev_seg)) {
        $cur_seg = $cur_seg->prev_seg;
        if (
            defined($cur_seg->low_end->bkpt) &&
            !$self->mate->eq($cur_seg->low_end->bkpt) &&
            !$cur_seg->low_end->bkpt->rg->is_small_td(%params) &&
            !$cur_seg->low_end->bkpt->rg->is_small_del(%params) &&
            !$cur_seg->low_end->bkpt->rg->is_part_of_shard_cycle(%params)
        ) {
            push @effector_rg_ends, $cur_seg->low_end->bkpt;
        }
    }

    return @effector_rg_ends;
}

sub get_ds_effector_rg_ends {
    my $self = shift;
    my %params = @_;
    my @effector_rg_ends = ();

    if (
        $self->is_rev &&
        defined($self->segment->high_end->bkpt) &&
        !$self->mate->eq($self->segment->high_end->bkpt)
    ) {
        push @effector_rg_ends, $self->segment->high_end->bkpt;
    }

    my $cur_seg = $self->segment;
    while (defined($cur_seg->next_seg)) {
        $cur_seg = $cur_seg->next_seg;
        if (
            defined($cur_seg->high_end->bkpt) &&
            !$self->mate->eq($cur_seg->high_end->bkpt) &&
            !$cur_seg->high_end->bkpt->rg->is_small_td(%params) &&
            !$cur_seg->high_end->bkpt->rg->is_small_del(%params) &&
            !$cur_seg->high_end->bkpt->rg->is_part_of_shard_cycle(%params)
        ) {
            push @effector_rg_ends, $cur_seg->high_end->bkpt;
        }
    }

    return @effector_rg_ends;
}

sub neighbours_array {
    # This subroutine looks for shard neighbours
    my $self = shift;
    my %params = @_;
    my @neighbours = ();
    if ($self->is_fwd) {
        if (
            defined($self->segment->next_seg) &&
            !defined($self->segment->next_seg->low_end->bkpt) &&
            $self->segment->next_seg->length <= $params{shard_bypassing_slop} &&
            defined($self->segment->next_seg->high_end->bkpt)
        ) {
            push @neighbours, $self->segment->next_seg->high_end->bkpt;
        }

        if (!$self->segment->is_bal_rg_gap(%params)) {
            if (
                defined($self->segment->prev_seg) &&
                !defined($self->segment->low_end->bkpt) &&
                $self->segment->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->prev_seg->high_end->bkpt)
            ) {
                push @neighbours, $self->segment->prev_seg->high_end->bkpt;
            }
        }
        else {
            # Difference is that now segment can have a low end breakpoint
            if (
                defined($self->segment->prev_seg) &&
                $self->segment->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->prev_seg->high_end->bkpt)
            ) {
                push @neighbours, $self->segment->prev_seg->high_end->bkpt;
            }
        }
    }
    else {
        if (
            defined($self->segment->prev_seg) &&
            !defined($self->segment->prev_seg->high_end->bkpt) &&
            $self->segment->prev_seg->length <= $params{shard_bypassing_slop} &&
            defined($self->segment->prev_seg->low_end->bkpt)
        ) {
            push @neighbours, $self->segment->prev_seg->low_end->bkpt;
        }

        if (!$self->segment->is_bal_rg_gap(%params)) {
            if (
                defined($self->segment->next_seg) &&
                !defined($self->segment->high_end->bkpt) &&
                $self->segment->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->next_seg->low_end->bkpt)
            ) {
                push @neighbours, $self->segment->next_seg->low_end->bkpt;
            }
        }
        else {
            # Difference is that now segment can have a low end breakpoint
            if (
                defined($self->segment->next_seg) &&
                $self->segment->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->next_seg->low_end->bkpt)
            ) {
                push @neighbours, $self->segment->next_seg->low_end->bkpt;
            }
        }
    }
    return @neighbours;
}

sub closest_us_rg_end {
    my $self = shift;
    my %params = @_;
    if (!exists($params{'within'})) {
        die "Parameter 'within' is required for $self\->closest_us_rg_end()";
    }
    if (!exists($params{'exclude_rgs'})) {
        die "Parameter 'exclude_rgs' is required for $self\->closest_us_rg_end()";
        # $params{exclude_rgs} = [];
    }

    if (
        $self->is_fwd &&
        defined($self->segment->low_end->bkpt) &&
        !(grep { $self->segment->low_end->bkpt->rg->id eq $_ } @{$params{exclude_rgs}})
    ) {
        if ($self->segment_end->pos - $self->segment_end->segment->low_end->pos > $params{'within'}) {
            return undef;
        }
        return $self->segment_end->segment->low_end->bkpt;
    }

    my $cur_seg = $self->segment_end->segment;
    my $rg_end;
    while (1) {
        last if !defined($cur_seg->prev_seg);
        $cur_seg = $cur_seg->prev_seg;

        if (defined($cur_seg->high_end->bkpt)) {
            $rg_end = $cur_seg->high_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($self->segment_end->pos - $rg_end->segment_end->pos > $params{'within'}) {
                    return undef;
                }
                return $rg_end;
            }
        }
        if (defined($cur_seg->low_end->bkpt)) {
            $rg_end = $cur_seg->low_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($self->segment_end->pos - $rg_end->segment_end->pos > $params{'within'}) {
                    return undef;
                }
                return $rg_end;
            }
        }
    }

    return undef;
}

sub closest_ds_rg_end {
    my $self = shift;
    my %params = @_;
    if (!exists($params{'within'})) {
        die "Parameter 'within' is required for $self\->closest_ds_rg_end()";
    }
    if (!exists($params{'exclude_rgs'})) {
        die "Parameter 'exclude_rgs' is required for $self\->closest_ds_rg_end()";
        # $params{exclude_rgs} = [];
    }

    if (
        $self->is_rev &&
        defined($self->segment->high_end->bkpt) &&
        !(grep { $self->segment->high_end->bkpt->rg->id eq $_ } @{$params{exclude_rgs}})
    ) {
        if ($self->segment_end->segment->high_end->pos - $self->segment_end->pos > $params{'within'}) {
            return undef;
        }
        return $self->segment_end->segment->high_end->bkpt;
    }

    my $cur_seg = $self->segment_end->segment;
    my $rg_end;
    while (1) {
        last if !defined($cur_seg->next_seg);
        $cur_seg = $cur_seg->next_seg;

        if (defined($cur_seg->low_end->bkpt)) {
            $rg_end = $cur_seg->low_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($rg_end->segment_end->pos - $self->segment_end->pos > $params{'within'}) {
                    return undef;
                }
                return $rg_end;
            }
        }
        if (defined($cur_seg->high_end->bkpt)) {
            $rg_end = $cur_seg->high_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($rg_end->segment_end->pos - $self->segment_end->pos > $params{'within'}) {
                    return undef;
                }
                return $rg_end;
            }
        }
    }

    return undef;
}


sub closest_rg_end {
    my $self = shift;
    my %params = @_;
    my $closest_us_rg_end = $self->closest_us_rg_end(%params);
    my $closest_ds_rg_end = $self->closest_ds_rg_end(%params);
    if (!defined($closest_us_rg_end) && !defined($closest_ds_rg_end)) {
        return undef;
    }
    elsif (!defined($closest_us_rg_end)) {
        return $closest_ds_rg_end;
    }
    elsif (!defined($closest_ds_rg_end)) {
        return $closest_us_rg_end;
    }
    else {
        if (
            $self->segment_end->pos - $closest_us_rg_end->segment_end->pos <=
            $closest_ds_rg_end->segment_end->pos - $self->segment_end->pos
        ) {
            return $closest_us_rg_end;
        }
        else {
            return $closest_ds_rg_end;
        }
    }
}
#
# End of worker subroutines
#


1;
