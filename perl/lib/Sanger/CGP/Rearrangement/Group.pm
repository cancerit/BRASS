package Sanger::CGP::Rearrangement::Group;

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

use warnings FATAL => 'all';
use strict;
use Scalar::Util qw(blessed);
use Math::Combinatorics;
use Bio::Brass qw($VERSION);

sub new {
    my $class = shift;
    my %params = @_;

    my $rg_group = bless {
        path => (exists($params{path}) ? $params{path} : undef),
        target_component => (
            exists($params{target_component}) ?
            $params{target_component} :
            undef
        ),
    }, $class;

    if (!exists($params{components})) {
        die "Argument 'components' is required for $class\->new()";
    }
    for (@{$params{components}}) {
        if (!defined(Scalar::Util::blessed $_) || Scalar::Util::blessed($_) ne 'Sanger::CGP::Rearrangement::ConnectedComponent') {
            die "The 'components' argument in $class\->new() must be of type 'Sanger::CGP::Rearrangement::ConnectedComponent'";
        }
        $rg_group->{components}->{$_} = $_;
    }

    return $rg_group;
}

#
# Getter methods
#
sub path {
    my $self = shift;
    return $self->{path};
}

sub components {
    my $self = shift;
    return $self->{components};
}

sub target_component {
    my $self = shift;
    if (!defined($self->{target_component})) {
        die "Attempted to access undefined $self\->{target_component} throught $self\->target_component()";
    }
    return $self->{target_component};
}

sub components_array {
    my $self = shift;
    return(values %{$self->components});
}

sub rgs {
    my $self = shift;
    my %rgs = ();
    for my $component ($self->components_array) {
        for (keys %{$component->rgs}) {
            $rgs{$_} = $component->rgs->{$_};
        }
    }
    return \%rgs;
}

sub rgs_array {
    my $self = shift;
    return(values %{$self->rgs});
}

sub size {
    my $self = shift;
    return scalar($self->components_array);
}

sub print {
    my $self = shift;
    print "$self, number of components: " . $self->size  . "\n";
    if (defined($self->path)) {
        $self->path->print;
    }
    else {
        print "Path: none\n";
    }
    for ($self->components_array) {
        $_->print;
    }
}

sub component_by_idx {
    my $self = shift;
    my $idx = shift;
    if (!defined($idx)) {
        return ($self->components_array)[0];
    }
    if ($idx < 0 or $idx > scalar($self->components_array) - 1) {
        die;
    }
    return ($self->components_array)[$idx];
}

sub chrs {
    my $self = shift;
    my %chrs = ();
    for my $component ($self->components_array) {
        %chrs = (%chrs, %{$component->chrs});
    }
    return \%chrs;
}

sub chrs_array {
    my $self = shift;
    return(values %{$self->chrs});
}

sub is_intra_chromosomal {
    my $self = shift;
    my @chrs = $self->chrs_array;
    for (@chrs[1..$#chrs]) {
        return 0 if $chrs[0] ne $_;
    }
    return 1;
}
#
# End of getter methods
#


#
# Worked methods
#
sub get_normalised_rg_patterns {
    my $self = shift;
    my %params = @_;
    my @involved_chrs = $self->chrs_array;
    my $involved_rgs = $self->rgs;
    my($cur_chr, $seg, %involved_segs_of_chr);

    for $cur_chr (@involved_chrs) {
        die if !$cur_chr->is_sorted;

        for $seg ($cur_chr->segments_array) {
            if (
                (   defined($seg->low_end->bkpt) &&
                    exists($involved_rgs->{$seg->low_end->bkpt->id}) ) ||
                (   defined($seg->prev_seg) &&
                    defined($seg->prev_seg->high_end->bkpt) &&
                    exists($involved_rgs->{$seg->prev_seg->high_end->bkpt->id}) ) ||
                (   defined($seg->high_end->bkpt) &&
                    exists($involved_rgs->{$seg->high_end->bkpt->id}) ) ||
                (   defined($seg->next_seg) &&
                    defined($seg->next_seg->low_end->bkpt) &&
                    exists($involved_rgs->{$seg->next_seg->low_end->bkpt->id}) )
            ) {
                push @{$involved_segs_of_chr{$cur_chr->name}}, $seg;
            }
        }
    }

    ## Rearrangements on all affected chromosomes can be represented on both
    ## orientations.
    my @permutations = Math::Combinatorics::permute(0..$#involved_chrs);
    my @chr_is_inverted = ();  ## A binary vector indicating whether a
                               ## chromosome is inverted in the current
                               ## permutation.
    my @cur_chr_ordering = ();
    my(@idx, $i, $j, $seg_idx, @segs_per_chr, $cur_rg_end, $facing_rg_end);
    ## Rearrangements are stored in 'buckets' sorted by when they are first
    ## encountered in each chromosomal arrangement.
    my(@rg_buckets, %bucket_idx_of_rg, @cn_change_across_bkpt, @cn_change_relative_var);
    my($rg_string, $segs_per_chr_string);
    my @best_patterns = ();

    ## For each order of chromosomes...
    for (0..$#permutations) {
        @cur_chr_ordering = @involved_chrs[@{$permutations[$_]}];

        ## And for each combination of orders...
        for (0..(2**scalar(@involved_chrs)-1)) {
            @chr_is_inverted = split '', sprintf("%0" . scalar(@cur_chr_ordering) . 'b', $_);
            @segs_per_chr = ();  ## Counts number of segments in each chromosome in current ordering.
            @rg_buckets = %bucket_idx_of_rg = ();
            $seg_idx = 0;  ## Index for every segment (across all chromosomes)
            @cn_change_across_bkpt = @cn_change_relative_var = ();

            ## And for each permuted chromosome order and orientation...
            for $i (0..$#cur_chr_ordering) {
                ## Some initiations
                push @segs_per_chr, 0;
                $cur_chr = $cur_chr_ordering[$i];
                if (!$chr_is_inverted[$i]) {
                    $j = 0;
                    while ($j <= $#{$involved_segs_of_chr{$cur_chr->name}}) {
                        $seg = $involved_segs_of_chr{$cur_chr->name}->[$j];

                        # Deal with low end
                        if (defined($seg->low_end->bkpt) && exists($involved_rgs->{$seg->low_end->bkpt->id})) {
                            $cur_rg_end = $seg->low_end->bkpt;
                            if ($cur_rg_end->end eq 'high' && $cur_rg_end->is_foldback(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }

                            if (exists($bucket_idx_of_rg{$cur_rg_end->id})) {
                                if (scalar(@{$rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]}) != 1) {
                                    die;
                                }
                                $rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]->[1] = "$seg_idx-";
                            }
                            else {
                                $rg_buckets[$#rg_buckets+1] = ["$seg_idx-"];
                                $bucket_idx_of_rg{$cur_rg_end->id} = $#rg_buckets;
                            }

                            $_ = $cur_rg_end->rg->weighted_avg_cn_change_across_rg(%params);
                            push @cn_change_across_bkpt, (defined($_) ? $_ : undef);
                            push @cn_change_relative_var, (defined($_) ? $cur_rg_end->rg->weighted_avg_cn_change_var_across_rg(%params) : 0);
                        }

                        # Then deal with high end in the same way
                        if (defined($seg->high_end->bkpt) && exists($involved_rgs->{$seg->high_end->bkpt->id})) {
                            $cur_rg_end = $seg->high_end->bkpt;
                            if ($seg->is_bal_rg_overlap(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }
                            if ($cur_rg_end->end eq 'high' && $cur_rg_end->is_foldback(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }

                            if (exists($bucket_idx_of_rg{$cur_rg_end->id})) {
                                if (scalar(@{$rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]}) != 1) {
                                    die;
                                }
                                $rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]->[1] = "$seg_idx+";
                            }
                            else {
                                $rg_buckets[$#rg_buckets+1] = ["$seg_idx+"];
                                $bucket_idx_of_rg{$cur_rg_end->id} = $#rg_buckets;
                            }

                            $_ = $cur_rg_end->rg->weighted_avg_cn_change_across_rg(%params);
                            push @cn_change_across_bkpt, (defined($_) ? -$_ : undef);
                            push @cn_change_relative_var, (defined($_) ? $cur_rg_end->rg->weighted_avg_cn_change_var_across_rg(%params) : 0);
                        }

                        if (
                            (
                                (
                                    defined($seg->high_end->bkpt) &&
                                    exists($involved_rgs->{$seg->high_end->bkpt->id})
                                ) ||
                                (
                                    $j+1 <= $#{$involved_segs_of_chr{$cur_chr->name}} &&
                                    defined($involved_segs_of_chr{$cur_chr->name}->[$j+1]->low_end->bkpt) &&
                                    exists($involved_rgs->{$involved_segs_of_chr{$cur_chr->name}->[$j+1]->low_end->bkpt->id})
                                )
                            ) &&
                            !$seg->is_bal_rg_gap(%params)
                        ) {
                            $seg_idx++;
                            $segs_per_chr[-1]++;
                        }
                        $j++;
                    }
                }
                else {
                    $j = $#{$involved_segs_of_chr{$cur_chr->name}};
                    while ($j >= 0) {
                        $seg = $involved_segs_of_chr{$cur_chr->name}->[$j];

                        # Deal with high end
                        if (defined($seg->high_end->bkpt) && exists($involved_rgs->{$seg->high_end->bkpt->id})) {
                            $cur_rg_end = $seg->high_end->bkpt;
                            if ($cur_rg_end->end eq 'low' && $cur_rg_end->is_foldback(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }

                            if (exists($bucket_idx_of_rg{$cur_rg_end->id})) {
                                if (scalar(@{$rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]}) != 1) {
                                    die;
                                }
                                $rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]->[1] = "$seg_idx-";
                            }
                            else {
                                $rg_buckets[$#rg_buckets+1] = ["$seg_idx-"];
                                $bucket_idx_of_rg{$cur_rg_end->id} = $#rg_buckets;
                            }

                            $_ = $cur_rg_end->rg->weighted_avg_cn_change_across_rg(%params);
                            push @cn_change_across_bkpt, (defined($_) ? $_ : undef);
                            push @cn_change_relative_var, $cur_rg_end->rg->weighted_avg_cn_change_var_across_rg(%params);
                        }

                        # Then deal with low end in the same way
                        if (defined($seg->low_end->bkpt) && exists($involved_rgs->{$seg->low_end->bkpt->id})) {
                            $cur_rg_end = $seg->low_end->bkpt;
                            if ($seg->is_bal_rg_overlap(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }
                            if ($cur_rg_end->end eq 'low' && $cur_rg_end->is_foldback(%params)) {
                                $seg_idx--;
                                $segs_per_chr[-1]--;
                            }

                            if (exists($bucket_idx_of_rg{$cur_rg_end->id})) {
                                if (scalar(@{$rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]}) != 1) {
                                    die;
                                }
                                $rg_buckets[$bucket_idx_of_rg{$cur_rg_end->id}]->[1] = "$seg_idx+";
                            }
                            else {
                                $rg_buckets[$#rg_buckets+1] = ["$seg_idx+"];
                                $bucket_idx_of_rg{$cur_rg_end->id} = $#rg_buckets;
                            }

                            $_ = $cur_rg_end->rg->weighted_avg_cn_change_across_rg(%params);
                            push @cn_change_across_bkpt, (defined($_) ? -$_ : undef);
                            push @cn_change_relative_var, $cur_rg_end->rg->weighted_avg_cn_change_var_across_rg(%params);
                        }

                        if (
                            (
                                (
                                    defined($seg->low_end->bkpt) &&
                                    exists($involved_rgs->{$seg->low_end->bkpt->id})
                                ) ||
                                (
                                    $j-1 >= 0 &&
                                    defined($involved_segs_of_chr{$cur_chr->name}->[$j-1]->high_end->bkpt) &&
                                    exists($involved_rgs->{$involved_segs_of_chr{$cur_chr->name}->[$j-1]->high_end->bkpt->id})
                                )
                            ) &&
                            !$seg->is_bal_rg_gap(%params)
                        ) {
                            $seg_idx++;
                            $segs_per_chr[-1]++;
                        }
                        $j--;
                    }
                }

                $seg_idx++;
                $segs_per_chr[-1]++;
            }

            ## After going through all segments in the current chromosome
            ## ordering/orientation, make the rearrangement string and
            ## see whether it's lexicographically smaller than the previous
            ## smallest rearrangement strings.
            @rg_buckets = sort { substr($a->[0], 0, -1) <=> substr($b->[0], 0, -1) || substr($b->[0], -1) cmp substr($a->[0], -1) } @rg_buckets;
            $rg_string = join(
                '/',
                map( {$_->[0] . ',' . $_->[1]} @rg_buckets )
            );
            $segs_per_chr_string = join('/', @segs_per_chr);

            ## If the current rearrangement string is lexicographically
            ## smaller than all the earlier ones then replace
            ## @best_patterns.
            if (
                @best_patterns == 0 ||
                (
                    "$rg_string $segs_per_chr_string" lt
                    ($best_patterns[0]->{rg_string} . ' ' . $best_patterns[0]->{segs_per_chr_string})
                )
            ) {
                @best_patterns = ({
                    rg_string => $rg_string,
                    segs_per_chr_string => $segs_per_chr_string,
                    chr_ordering => [@cur_chr_ordering],
                    chr_is_inverted => [@chr_is_inverted],
                    cn_changes => [@cn_change_across_bkpt],
                    cn_change_variances => [@cn_change_relative_var],
                });
            }
            elsif (
                "$rg_string $segs_per_chr_string" eq
                ($best_patterns[0]->{rg_string} . ' ' . $best_patterns[0]->{segs_per_chr_string})
            ) {
                push @best_patterns, {
                    rg_string => $rg_string,
                    segs_per_chr_string => $segs_per_chr_string,
                    chr_ordering => [@cur_chr_ordering],
                    chr_is_inverted => [@chr_is_inverted],
                    cn_changes => [@cn_change_across_bkpt],
                    cn_change_variances => [@cn_change_relative_var],
                };
            }
        }
    }

    $self->{best_patterns} = [@best_patterns];
}

sub print_normalised_rg_patterns {
    my $self = shift;
    if (!exists($self->{best_patterns}) || !defined($self->{best_patterns})) {
        die "Attempted to call $self\->print_normalised_rg_patterns() with undefined $self\->{best_patterns}";
    }

    for my $p (@{$self->{best_patterns}}) {
        print('rg pattern: ' . $p->{rg_string} . ' ' . $p->{segs_per_chr_string} . " (rearrangements segments_per_chr)\n");
        print '  ' . join(' ', map( { if (!defined) { ' NA' } else { sprintf('% .3f', $_) } } @{$p->{cn_changes}})) . " (cn change over rg)\n";
        print '  ' . join(' ', map( { if (!defined) { ' NA' } else { sprintf('% .3f', $_) } } @{$p->{cn_change_variances}})) . " (variance of cn change)\n";
    }
}

sub unique_string {
    my $self = shift;
    return join(
        '+',
        $self->target_component,
        sort(grep {$_ ne $self->target_component} keys(%{$self->components})),
    );
}
#
# End of worker methods
#


1;
