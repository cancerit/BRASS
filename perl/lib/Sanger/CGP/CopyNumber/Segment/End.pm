package Sanger::CGP::CopyNumber::Segment::End;

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
