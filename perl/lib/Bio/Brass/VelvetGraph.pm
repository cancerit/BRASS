package Bio::Brass::VelvetGraph;

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


use strict;
use warnings;

use base qw(Exporter);
our @EXPORT = qw(add_graph_labels set_vertex_style view_graph map_contigs);

use File::Temp qw(tempfile);

=head1 NAME

Bio::Brass::VelvetGraph - Velvet graph manipulation

=head1 SYNOPSIS

    use Graph;
    use Bio::Brass::VelvetGraph;

    add_graph_labels($graph);
    add_graph_labels($graph, $vlabel, $vtooltip);

    set_vertex_style($graph, "filled", @quintet);

    view_graph($graph, $filename);

    map_contigs($exonerate_opts, $graph, @vertices);

=head1 DESCRIPTION

=cut

use Graph;
use Graph::Reader::Velvet; # for get_vertex_contig()
use Graph::Writer::Dot;
use Bio::Brass::Alignment;
use Bio::Brass qw($VERSION);

=head2 add_graph_labels

Adds attributes for consumption by B<neato>.  The C<$vlabel> and C<$vtooltip>
template arguments are optional, and may contain the following conversion
specifications: C<%I> (id), C<%L> (contig length), C<%C> (contig).
If either or both of these template arguments is/are missing,
the corresponding vertex attributes will not be set.

=cut

sub add_graph_labels {
    my ($graph, $label, $tooltip) = @_;

    $graph->set_graph_attribute('layout', 'neato');
    #$graph->set_graph_attribute("name", "whatever"); # FIXME Not listened to

    my ($need_length, $need_contig);

    foreach ($label, $tooltip) {
	next unless defined $_;
	s'%I'%1$d'g;
	s'%L'%2$d'g and $need_length = 1;
	s'%C'%3$s'g and $need_contig = 1;
    }

    foreach ($graph->vertices) {
	my $len = $need_length? $graph->get_vertex_contig_length($_) : 0;
	my $ctg = $need_contig? $graph->get_vertex_contig($_) : '?';

	$graph->set_vertex_attribute($_, 'label',
	    sprintf($label, $_, $len, $ctg)) if defined $label;

	$graph->set_vertex_attribute($_, 'tooltip',
	    sprintf($tooltip, $_, $len, $ctg)) if defined $tooltip;
    }

    foreach ($graph->edges) {
	$graph->set_edge_attribute(@$_, 'label', $graph->get_edge_weight(@$_));
	$graph->set_edge_attribute(@$_, 'fontsize', 10);
    }
}

=head2 set_vertex_style

Sets the C<style> attribute as specified for each of the specified vertices.

=cut

sub set_vertex_style {
    my $graph = shift @_;
    my $style = shift @_;

    $graph->set_vertex_attribute($_, 'style', $style) foreach @_;
}

=head2 view_graph

Displays the graph on screen, as laid out by B<neato>.
Invokes C<add_graph_labels($graph, "%I %L", "%C")> if this has not
already been called for this graph.

=cut

sub view_graph {
    my ($graph, $fname) = @_;

    die 'FIXME: temporary view_graph filename' unless defined $fname;

    add_graph_labels($graph, '%I %L', '%C')
	    unless $graph->has_graph_attribute('layout');

    my $writer = Graph::Writer::Dot->new();
    if ($^O =~ /darwin/i) {
	$writer->write_graph($graph, $fname);
	system("open -a Graphviz $fname") == 0 or die "Can't spawn open: $!\n";
	}
    else {
	open my $NEATO, '|neato -Tgtk' or die "Can't pipe to neato: $!\n";
	$writer->write_graph($graph, $NEATO);
	close $NEATO or die "neato execution failed\n";
    }
}

=head2 map_contigs

Invokes exonerate on the specified vertices' contigs, adding a C<mappings>
attribute to each vertex, in the form C<blah>.  Unless C<$options> specifies
a particular model, exonerate is invoked again on any unmappable vertices
using the I<affine:local> model, which allows for indels.

=cut

sub map_contigs {
    my ($working_dir, $options, $graph, $regions, @vertices) = @_;
    my @contigs = map { $graph->get_vertex_contig($_) } @vertices;

    # FIXME  Only really need to do pos/neg vertex pairs once.  Maybe.
    ## Assumes the mappings come out in the same order as the contigs went in...
    foreach my $mappings (map_sequences($working_dir, $options, $regions, @contigs)) {

        my $vertex = shift @vertices;
        my $contig = shift @contigs;
        @$mappings = sort @$mappings;
        $graph->set_vertex_attribute($vertex, 'mappings', $mappings);

        my $contig_length = length $contig;

        ## If no mappings found
        ## or the first mapping does not start at pos 1
        ## or the first mapping does not end at er... the end
        if (scalar(@$mappings) == 0 || $$mappings[0]->query->start != 1 ||
            $$mappings[0]->query->end != $contig_length) {
            push @vertices, $vertex;
            push @contigs, $contig;
        }
    }

    if ($options !~ /--model/ && scalar(@vertices) > 0) {
    	warn join(' ','Aligning with --model affine:local for vertices:',map{$_}@vertices);
        foreach my $mappings (map_sequences($working_dir, "$options --model affine:local", $regions, @contigs)) {
            my $vertex = shift @vertices;
            @$mappings = sort @$mappings;
            $graph->set_vertex_attribute($vertex, 'mappings', $mappings);
        }
    }
}

=head1 AUTHOR

John Marshall E<lt>jm18@sanger.ac.ukE<gt>

=cut

1;
