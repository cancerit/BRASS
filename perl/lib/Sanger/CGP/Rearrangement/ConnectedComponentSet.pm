package Sanger::CGP::Rearrangement::ConnectedComponentSet;

use warnings FATAL => 'all';
use strict;
use Bio::Brass qw($VERSION);

use Scalar::Util qw(blessed);
use Sanger::CGP::Rearrangement::ConnectedComponent;

sub new {
    my $class = shift;
    my %params = @_;
    my $rcc_set = bless { components => {} }, $class;

    if (exists $params{components}) {
        for (@{$params{components}}) {
            if (!defined(Scalar::Util::blessed $_) or Scalar::Util::blessed($_) ne 'Sanger::CGP::Rearrangement::ConnectedComponent') {
                die "Array elements of dereferenced argument 'components' must be of class 'Sanger::CGP::Rearrangement::ConnectedComponent' in $class\->Sanger::CGP::Rearrangement::ConnectedComponentSet";
            }
            if (!exists($rcc_set->{components}->{$_})) {
                $rcc_set->{components}->{$_} = $_;
            }
        }
    }

    return $rcc_set;
}

#
# Getter and setter methods
#
sub components {
    my $self = shift;
    return $self->{components};
}

sub components_array {
    my $self = shift;
    return values(%{$self->components});
}

sub print {
    my $self = shift;
    print "$self, number of components: " . $self->size . "\n";
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

sub add {
    my $self = shift;
    for (@_) {
        if (!defined(Scalar::Util::blessed $_) or Scalar::Util::blessed($_) ne 'Sanger::CGP::Rearrangement::ConnectedComponent') {
            die "Array elements of dereferenced argument 'components' must be of class 'Sanger::CGP::Rearrangement::ConnectedComponent' in $self\->Sanger::CGP::Rearrangement::ConnectedComponentSet";
        }
        if (!exists($self->components->{$_})) {
            $self->components->{$_} = $_;
        }
    }
}

sub delete {
    my $self = shift;
    my $target = shift;
    if (!exists($self->components->{$target})) {
        die "Attempted to remove nonexistent component '$target' in $self\->delete()";
    }
    delete($self->{components}->{$target});
}

sub first {
    my $self = shift;
    my ($first_key) = (sort { $a cmp $b } keys(%{$self->components}))[0];
    return $self->components->{$first_key};
}

sub delete_first {
    my $self = shift;
    my $first = $self->first;
    $self->delete($first);
}

sub shift_first {
    my $self = shift;
    my $first = $self->first;
    $self->delete($first);
    return $first;
}

sub copy {
    my $self = shift;
    return bless { components => {%{$self->components}} }, ref($self);
}

sub contains {
    my $self = shift;
    my $target = shift;
    if (!defined(Scalar::Util::blessed $target) or Scalar::Util::blessed($target) ne 'RearrangementConnectedComponent') {
        die "Array elements of dereferenced argument 'components' must be of class 'RearrangementConnectedComponent' in $self\->RearrangementConnectedComponentSet";
    }

    return exists($self->components->{$target});
}

sub size {
    my $self = shift;
    return scalar(keys %{$self->components});
}
#
# End of getter and setter methods
#

1;
