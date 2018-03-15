package Graph::Reader::Velvet;

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

use base qw(Graph::Reader);
our $VERSION = '0.83';

=head1 NAME

Graph::Reader::Velvet - read a Velvet graph file into a Graph instance

=head1 SYNOPSIS

    use Graph;
    use Graph::Reader::Velvet;

    $reader = Graph::Reader::Velvet->new();
    $graph = $reader->read_graph('LastGraph');

    $k = $graph->get_graph_hash_length();
    $seq = $graph->get_vertex_seq($vertex);
    $contig = $graph->get_vertex_contig($vertex);
    $length = $graph->get_vertex_contig_length($vertex);

=cut

use Graph;

sub _read_graph {
    my ($self, $G, $FILE) = @_;

    my (undef, $count, $k) = split /\s+/, <$FILE>;
    $G->set_graph_attribute('hash_length', $k);
    $G->set_graph_attribute('sequence_count', $count);

    my $unget = undef;
    while ($_ = (defined $unget)? $unget : <$FILE>) {
		undef $unget;
		chomp $_;

		if (/^NODE\s/) {
	        my (undef, $id, $node_len_in_kmers, @coverage) = split /\s+/, $_;
		    $id = int($id);

		    my ($seq, $twinseq);
		    if (scalar @coverage > 0) {
			# An ordinary file, listing both sequences explicitly
			chomp ($seq = <$FILE>);
			chomp ($twinseq = <$FILE>);
		    }
		    else {
			# A PreGraph file, with the sequences conjoined
			chomp ($_ = <$FILE>);
			$seq = substr($_, $k-1);
			$twinseq = reverse substr($_, 0, -($k-1));
			$twinseq =~ tr/ACGTacgt/TGCAtgca/;
		    }

		    $G->add_vertex($id);
		    $G->add_vertex(-$id);
		    $G->set_vertex_attribute($id,  'seq', $seq);
		    $G->set_vertex_attribute(-$id, 'seq', $twinseq);

		    my $cov_total = 0;
		    for (my $cat = 1; scalar @coverage > 0; $cat++) {
			my $cov  = shift @coverage;
			my $ocov = shift @coverage;
			$G->set_vertex_attribute($id,  "short${cat}_cov",  $cov);
			$G->set_vertex_attribute($id,  "short${cat}_Ocov", $ocov);
			$G->set_vertex_attribute(-$id, "short${cat}_cov",  $cov);
			$G->set_vertex_attribute(-$id, "short${cat}_Ocov", $ocov);
			$cov_total += $cov;
		    }

		    $G->set_vertex_weight($id,  $cov_total);
		    $G->set_vertex_weight(-$id, $cov_total);
		}
		elsif (/^ARC\s/) {
		    my (undef, $start, $end, $multiplicity) = split /\s+/, $_;
		    $start = int($start);
		    $end = int($end);

		    $G->add_weighted_edge($start, $end, $multiplicity);
		    $G->add_weighted_edge(-$end, -$start, $multiplicity);
		}
		elsif (/^SEQ\s/) {
		    my (undef, $seqId) = split /\s+/, $_;
		    while (<$FILE>) {
				last if $_ !~ /^-?\d/;
				chomp;
				my ($id, $offset, $start, $end, $endOffset) = split /\s+/, $_;
				my $seq_pos = $G->get_vertex_attribute($id,  'seq_positions') || {};
				$G->set_vertex_attribute($id,  'seq_positions', $seq_pos) unless(scalar keys %$seq_pos);
				push(@{$seq_pos->{$seqId}},[$offset, $start, $end, $endOffset]);
		    }
		    $unget = $_;
		}
		elsif (/^NR\s/) {
		    my (undef, $id, $shortReadCount) = split /\s+/, $_;

		    $G->set_vertex_attribute($id,  'short_read_count', $shortReadCount);
			$G->set_vertex_attribute(-$id, 'short_read_count',  $shortReadCount);

		    my $short_reads = [];
		    $G->set_vertex_attribute($id, 'short_reads', $short_reads);
		    $G->set_vertex_attribute(-$id, 'short_reads', $short_reads);

		    while (<$FILE>) {
				last if $_ !~ /^-?\d/;
				chomp;
				#my ($readId, $offset, $start) = split /\s+/, $_;
				push(@$short_reads,[split /\s+/, $_]);
		    }
		    $unget = $_;
		}
		else {
		    warn qq{unknown velvet line "$_"\n};
		}
    }

    return 1;
}

=head1 DESCRIPTION

B<Graph::Reader::Velvet> is a class for reading in a directed graph in the
file format produced by I<Velvet>, a I<de novo> genomic assembler designed
for short read sequencing technologies.

The constructor, C<new()>, takes no arguments and generates a new reader
instance.  The C<read_graph()> method takes one argument, which may be either
a filename or a previously opened filehandle, and returns a directed graph.
B<Graph::Reader::Velvet> is a subclass of B<Graph::Reader>, which should be
consulted for further details of these generic Graph reader methods.

Attributes attached to the returned digraph represent its biological data.
The digraph itself has C<hash_length> and C<sequence_count> graph attributes,
containing the hash length (I<k> value) and sequence count fields from the
header line of the Velvet graph file.

The digraph has two vertices, I<N> and I<-N>, for each node/twin-node block.
Each vertex has a C<seq> attribute containing its sequence, i.e., the sequence
of the final nucleotides of its I<k>-mers.  Except for PreGraphs, each vertex
also has C<short1_cov>, C<short1_Ocov>, C<short2_cov>, and C<short2_Ocov>
attributes (and so on, if Velvet has been built with more categories)
containing its I<k>-mer coverage values, and is weighted with the total
of the C<shortN_cov> coverage values, i.e., it has a C<weight> attribute,
which can be retrieved with C<get_vertex_weight()>.

Each edge is weighted according to its multiplicity, i.e., it has a C<weight>
attribute, which can be retrieved with C<get_edge_weight()>.

=cut

{
package Graph;

=pod

For convenience, this package also adds methods to the Graph class
for retrieving the hash length and vertices' sequences and contigs.

=over 4

=item get_graph_hash_length

=item get_vertex_seq

These return the corresponding graph or vertex attribute.

=cut

sub get_graph_hash_length {
    my ($G) = @_;
    return $G->get_graph_attribute('hash_length');
}

sub get_vertex_seq {
    my $G = shift @_;
    return $G->get_vertex_attribute(@_, 'seq');
}

=item get_vertex_contig

Returns the vertex's contig, reconstructing it in the same way as Velvet
itself, by prepending a prefix of the reverse complement of the twin vertex's
sequence, and additionally, if the vertex's sequence is shorter than the
graph's hash length, inserting C<N>s S<in between> as appropriate.

=cut

sub get_vertex_contig {
    my ($G, $id) = @_;
    $id = int($id);
    my $k = $G->get_graph_hash_length();

    my $seq = $G->get_vertex_seq($id);
    return unless defined $seq;

    # The first $k bases of the reverse complement of the twin's sequence
    local $_ = reverse substr($G->get_vertex_seq(-$id), -$k);
    tr/ACGTacgt/TGCAtgca/;

    if (length $_ < $k) {
	my $n = /[ACGTN]/? 'N' : 'n';
	$_ .= $n x ($k - length $_);
    }

    my $overlap = lc substr($_, -1);
    warn "Badly overlapped vertex contig: $_|$seq\n"
	unless $overlap eq lc substr($seq, 0, 1) || $overlap eq 'n';

    return substr($_, 0, -1) . $seq;
}

=item get_vertex_contig_length

Returns the length of the vertex's contig, without actually reconstructing it.

=cut

sub get_vertex_contig_length {
    my $G = shift @_;
    my $seq = $G->get_vertex_seq(@_);
    return unless defined $seq;
    return length($seq) + $G->get_graph_hash_length() - 1;
}

} # End of package Graph additions

1;
__END__

=back

=head1 CAVEATS

Colorspace is not supported.  If the input graph file was produced
by S<a colorspace> version of Velvet, C<get_vertex_contig()> and
F<PreGraph> reading will be incorrect due to inappropriate
reverse complementing.

At present, B<Graph::Reader::Velvet> does nothing with C<SEQ> and C<NR>
stanzas in the input file.

Reading a F<PreGraph> returns an empty digraph with annotated vertices but
no edges, as PreGraph edges are regenerated by Velvet rather than being
stored in the file.

=head1 SEE ALSO

=over 4

=item L<http://www.ebi.ac.uk/~zerbino/velvet/>

The Velvet home page, with links to source code, documentation, papers,
and mailing list.

=item L<Graph>

Jarkko Hietaniemi's classes for representing graphs.

=item L<Graph::Reader>

The base class for B<Graph::Reader::Velvet>, containing definitions and
documentation of the generic C<new()> and C<read_graph()> reader methods.

=back

=head1 AUTHOR

John Marshall E<lt>john.marshall@sanger.ac.ukE<gt>

=cut

