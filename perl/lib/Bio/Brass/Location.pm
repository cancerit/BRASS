package Bio::Brass::Location;

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


use strict;
use warnings;
use Bio::Brass qw($VERSION);

=head1 NAME

Bio::Brass::Location - represent a location something something

=head1 SYNOPSIS

    use Bio::Brass::Location;

    $location = Bio::Brass::Location->new(@arguments);
    $location2 = $location->clone();

    $locn->start, $locn->end, $locn->length, $locn->strand, $locn->name;

    $location->trim($left, $right);
    $location->untrim();
    $location->revcomp();

=head1 DESCRIPTION

=cut

use constant { _POS => 0, _LIM => 1, _SEQ => 2, _OFFSET => 3, _NAME => 4 };

use overload
    '""' => \&to_string,
    'cmp' => \&compare;

=head2 new

    $location = Bio::Brass::Location->new
                    ($seqtext, $begin, $end, $strand, $name);

$name is optional.

=cut

sub new {
    my ($class, $seq, $begin, $end, $strand, $name)  = @_;

    my @fields;
    if ($strand eq '+') { @fields = (0, $end - $begin, $seq, $begin + 1) }
    else { @fields = ($end - $begin, 0, $seq, $end) }
    push @fields, $name if defined $name;

    return bless [@fields], $class;
}

sub clone {
    my ($self) = @_;
    return bless [@$self], ref($self);
}

sub to_string {
    my ($self) = @_;
    my $len = CORE::length $$self[_SEQ];
    my $name = exists $$self[_NAME]? "$$self[_NAME]: " : q{};
    my ($start, $end) = ($self->start, $self->end);
    return "{$name$start--$end $$self[_OFFSET]+[$$self[_POS],$$self[_LIM]) ${len}bp}";
}

sub compare {
    my ($a, $b, $reversed) = @_;

    my $diff = $a->start <=> $b->start;
    $diff = $a->end <=> $b->end if $diff == 0;
    return $reversed? -$diff : +$diff;
}

sub _isforward { my ($self) = @_; return $$self[_POS] >= 0 }

sub start  { my ($self) = @_; return $$self[_OFFSET] + abs($$self[_POS]) }
sub end    { my ($self) = @_; return $$self[_OFFSET] + abs($$self[_LIM] - 1) }
sub length { my ($self) = @_; return $$self[_LIM] - $$self[_POS] }
sub strand { my ($self) = @_; return $self->_isforward()? '+' : '-' }
sub name   { my ($self) = @_; return $$self[_NAME] }

sub seq {
    my ($self) = @_;
    return substr($$self[_SEQ], $$self[_POS], $$self[_LIM] - $$self[_POS]);
}

sub prettyname {
    my ($self, $prefix) = @_;

    local $_ = (defined $prefix)? $prefix : 'Chr.';
    $_ .= $$self[_NAME];
    $_ .= '-' unless $self->_isforward();
    return $_;
}

=head2 trim

    $location->trim($left, $right);

Trims the specified amounts from each end.

=cut

sub trim {
    my ($self, $left, $right) = @_;
    $$self[_POS] += $left;
    $$self[_LIM] -= $right;
    return $self;
}

=head2 untrim

    $location->untrim();

Undoes any trimming and other cookery, returning $location to its original
boundaries.

=cut

sub untrim {
    my ($self) = @_;
    if ($self->_isforward()) {
	$$self[_POS] = 0;
	$$self[_LIM] = CORE::length($$self[_SEQ]);
    }
    else {
	$$self[_POS] = -CORE::length($$self[_SEQ]);
	$$self[_LIM] = 0;
    }

    return $self;
}

=head2 revcomp

    $location->revcomp();

Flip the location to the other strand, swapping start and end coordinates
and reverse complementing the sequence.

=cut

sub revcomp {
    my ($self) = @_;

    $$self[_OFFSET] += $self->_isforward()? -1 : +1;
    ($$self[_POS], $$self[_LIM]) = (-$$self[_LIM], -$$self[_POS]);
    ($$self[_SEQ] = reverse $$self[_SEQ]) =~ tr/ACGTacgt/TGCAtgca/;
}

=head1 AUTHOR

John Marshall E<lt>jm18@sanger.ac.ukE<gt>

=cut

1;
