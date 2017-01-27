package Sanger::CGP::Rearrangement::ConnectedComponent;

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

use Scalar::Util qw(blessed);
use Bio::Brass qw($VERSION);

sub new {
    my $class = shift;
    my %params = @_;

    return bless {
        rgs => ($params{rgs} or {}),
    }, $class;
}

#
# Getter method
#
sub rgs {
    my $self = shift;
    if (!defined($self->{rgs})) {
        die "Attempted to call rgs() to get undefined $self\->{rgs}!";
    }
    return $self->{rgs};
}

sub rgs_array {
    my $self = shift;
    return (values %{$self->rgs});
}

sub rg_by_idx {
    my $self = shift;
    my $idx = shift;
    if (!defined($idx)) {
        return ($self->rgs_array)[0];
    }
    if ($idx < -1 or $idx > scalar($self->rgs_array) - 1) {
        die;
    }
    return ($self->rgs_array)[$idx];
}
#
# End of getter methods
#

#
# Helper methods
#
sub add_rgs {
    my $self = shift;
    my @rgs = shift;
    for (@rgs) {
        if (Scalar::Util::blessed($_) ne 'Sanger::CGP::Rearrangement') {
            die "Attempted to add a non-Rearrangement object in $self\->add_rgs()";
        }
        if (exists($self->{rgs}->{$_->id})) {
            warn 'Trying to add pre-existing rearrangement ' . $_->id . " into $self\->{rgs}.\n";
        }
        $self->{rgs}->{$_->id} = $_;
    }
}

sub size {
    my $self = shift;
    return scalar(keys %{$self->rgs});
}

sub chrs {
    my $self = shift;
    my %chrs = ();
    for ($self->rgs_array) {
        $chrs{$_->low_end->chr_name}  = $_->low_end->chr;
        $chrs{$_->high_end->chr_name} = $_->high_end->chr;
    }
    return(\%chrs);
}

sub to_s {
    my $self = shift;
    my $out_string = "$self, number of rgs: " . $self->size . "\n";
    for (values %{$self->rgs}) {
        $out_string .= $_->to_s . "\n";
    }
    return $out_string;
}

sub print {
    my $self = shift;
    print $self->to_s;
}
#
# End of helper methods
#

#
# Worker methods
#
sub two_rgs_form_inverted_duplication {
    my $rg1 = shift;
    my $rg2 = shift;
    my %params = @_;

    ## In this function $rg1 is the rearrangement that forms a
    ## fold-back rearrangement.

    if ($rg1->is_foldback(%params)) {
        if (
            $rg1->low_end->is_fwd &&
            $rg1->high_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->high_end->balanced_bkpt_partner_rg_end(%params) == $rg2->high_end &&
            $rg2->low_end->chr_name eq $rg2->high_end->chr_name &&
            $rg2->low_end->dir eq $rg2->high_end->dir
        ) {
            return 1;
        }
        if (
            $rg1->low_end->is_rev &&
            $rg1->low_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->low_end->balanced_bkpt_partner_rg_end(%params) == $rg2->low_end &&
            $rg2->low_end->chr_name eq $rg2->high_end->chr_name &&
            $rg2->low_end->dir eq $rg2->high_end->dir
        ) {
            return 1;
        }
    }
    return 0;
}

sub naive_classification {
    my $self = shift;
    my %params = @_;

    if ($self->size == 1) {
        my($rg) = $self->rgs_array;

        ## Reciprocal shards
        if ($rg->is_part_of_shard_cycle(%params)) {
            return 'shard_cycle';
        }

        if ($rg->low_end->chr_name eq $rg->high_end->chr_name) {
            if ($rg->is_foldback(%params)) {
                return 'fb';
            }
            if ($rg->low_end->dir eq $rg->high_end->dir) {
                return 'fb_like';  ## Fix definition
            }
            if ($rg->low_end->is_fwd && $rg->high_end->is_rev) {
                return 'del';
            }
            if ($rg->low_end->is_rev && $rg->high_end->is_fwd) {
                return 'td';
            }
            die;  ## Shouldn't get this far
        }
        else {
            return 'ut';
        }
    }
    elsif ($self->size == 2) {
        my ($rg1, $rg2) = $self->rgs_array;

        ## Direct inversion
        if (
            $rg1->low_end->chr_name eq $rg1->high_end->chr_name &&
            $rg1->low_end->dir eq $rg1->high_end->dir &&
            $rg1->low_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->low_end->balanced_bkpt_partner_rg_end(%params) == $rg2->low_end &&
            $rg1->high_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->high_end->balanced_bkpt_partner_rg_end(%params) == $rg2->high_end
        ) {
            return 'inv';
        }

        ## Balanced translocation
        if (
            $rg1->low_end->chr_name ne $rg1->high_end->chr_name &&
            $rg1->low_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->low_end->balanced_bkpt_partner_rg_end(%params) == $rg2->low_end &&
            $rg1->high_end->has_balanced_bkpt_partner_rg_end(%params) &&
            $rg1->high_end->balanced_bkpt_partner_rg_end(%params) == $rg2->high_end
        ) {
            return 'bt';
        }

        ## The two different inverted duplications
        if (
            two_rgs_form_inverted_duplication($rg1, $rg2, %params) ||
            two_rgs_form_inverted_duplication($rg2, $rg1, %params)
        ) {
            return 'id';
        }

        return 'complex';
    }
    else {
        ## Automatically complex if the component consists of more than
        ## 2 rearrangements.
        return 'complex';
    }
}
#
# End of worker methods
#


1;
