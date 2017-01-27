package Sanger::CGP::CopyNumber::Segment::End;

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

# Every copy number segment is modelled as two reciprocal
# copy number segment ends, which are linked to each other.

sub new {
    my $class = shift;
    my %params = @_;
    return bless {
        pos      => $params{pos},
        end      => ($params{end} or undef),  # 'low' or 'high'
        bkpt     => ($params{bkpt} or undef),     # Breakpoint is a rearrangement breakpoint associated with a segment.
        boundary => ($params{boundary} or undef), # Boundary is the rearrangement or copy number segmentation breakpoint
                                                  # demarcating the copy number segment.
        segment  => ($params{segment} or undef),
    }, $class;
}

#
# Getter methods
#
sub pos {
    my $self = shift;
    if (!defined($self->{pos})) {
        die "Attempted to call pos() to get undefined Sanger::CGP::CopyNumber::Segment::End->{pos}!";
    }
    return $self->{pos};
}

sub end {
    my $self = shift;
    if (!defined($self->{end})) {
        die "Attempted to call end() to get undefined Sanger::CGP::CopyNumber::Segment::End->{end}!";
    }
    return $self->{end};
}

sub bkpt {
    my $self = shift;
    return $self->{bkpt};
}

sub boundary {
    my $self = shift;
    if (!defined($self->{boundary})) {
        die "Attempted to call boundary() to get undefined Sanger::CGP::CopyNumber::Segment::End->{boundary} at segment " . $self->segment->segment_position;
    }
    return $self->{boundary};
}

sub segment {
    my $self = shift;
    if (!defined($self->{segment})) {
        die "Attempted to call segment() to get undefined Sanger::CGP::CopyNumber::Segment::End->{segment}!";
    }
    return $self->{segment};
}
#
# End of getter methods
#

sub bounded_by {
    my $self = shift;
    if (defined($self->bkpt)) {
        return('rg');
    }
    elsif ($self->boundary =~ /tel$/) {
        return('telomere');
    }
    else {
        return($self->boundary);
    }
}
sub is_bounded_by_rg {
    my $self = shift;
    my $rg = shift;

    if (defined($rg)) {
        if ($self->bounded_by ne 'rg') {
            return 0;
        }
        elsif (ref($rg) eq 'SCALAR') {
            return $self->bkpt->id eq $rg;
        }
        elsif (defined(Scalar::Util::blessed($rg))) {
            ## Treat as Rearrangement or RearrangementEnd
            if (Scalar::Util::blessed($rg) ne 'Sanger::CGP::Rearrangement' && Scalar::Util::blessed($rg) ne 'Sanger::CGP::Rearrangement::End') {
                die;
            }
            return $self->bkpt->id == $rg->id;
        }
        else {
            return 0;
        }
    }

    return $self->bounded_by eq 'rg';
}
sub is_bounded_by_cn_bkpt { return $_[0]->bounded_by eq 'cn_bkpt'; }
sub is_bounded_by_tel     { return $_[0]->bounded_by eq 'telomere'; }
sub is_bounded_by_cen     { return $_[0]->bounded_by eq 'centromere'; }


1;
