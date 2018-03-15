package Bio::Brass;

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

use base qw(Exporter);
our @EXPORT = qw(find_breakpoints find_dusty_vertices
		 concatenate_path try_smooth_bubble update_bubble_consensus
		 smooth_fringe_bubbles annotate_bubbles
		 quintet_weights quintet_score
		 format_quintet format_quintet_alignment get_quintet_triplet
		 get_bp_region get_bp_surrounding_region
		 get_bp_alignment get_bp_alignment_text
		 is_dusty get_isolated_bp_alignment get_isolated_bp_surrounding_region
		 $VERSION);

our $VERSION = '6.1.0';

=head1 NAME

Bio::Brass - routines for identifying breakpoints

=head1 SYNOPSIS

    use Graph;
    use Bio::Brass;

    $dusty = is_dusty($graph);
    @quintets = find_breakpoints($graph);
    @vertices = find_dusty_vertices($graph);

    binmode STDOUT, ":utf8";
    print format_quintet(@quintet), "\n";
    print format_quintet_alignment($graph, @quintet);

    ($text, $width) = get_bp_region($graph, @quintet);
    $text = get_bp_alignment_text($graph, @quintet);
    @fields = get_bp_alignment($graph, @quintet);

=head1 DESCRIPTION

=cut

use Carp;
use Graph;
use Graph::Reader::Velvet; # for get_vertex_contig() et al
use List::Util qw(min max);

# Given a list of integers, returns how many of them are less than zero.
sub count_negatives {
    my $n = 0;
    foreach (@_) { $n++ if $_ < 0 }
    return $n;
}

=head2 is_dusty

    $dusty = is_dusty($graph);

Returns whether the graph consists sol

=cut

sub is_dusty {
    return 0;
}

=head2 find_dusty_vertices

    @vertices = find_dusty_vertices($graph);

Returns a list of the graph's isolated vertices that have
non-trivially long contigs.

=cut

sub find_dusty_vertices {
    my ($g, $factor) = @_;
    $factor = 2 unless defined $factor;

    my $min = $factor * $g->get_graph_hash_length();

    # FIXME or is_self_loop_vertex?
    return grep { $_ > 0 && $g->is_isolated_vertex($_) &&
		  $g->get_vertex_contig_length($_) >= $min } $g->vertices;
}

=head2 find_breakpoints

Single-copy breakpoints appear in Velvet's de Bruijn graph as a 5-vertex path
in the pattern S<C<v0 E<lt>- v1 -E<gt> v2 -E<gt> v3 E<lt>- v4>>.

    @quintets = find_breakpoints($graph);

Returns a list of references to 5-element arrays, each containing an instance
of the pattern.

=cut

sub find_breakpoints {
    my ($g) = @_;
    local $_;

    my $k = $g->get_graph_hash_length();

    # FIXME Previously had ">= 2", but how to decide which of the
    # odd ones out are the path endpoints?
    my @two_in  = grep { $g->in_degree($_) == 2 &&
			 length $g->get_vertex_seq($_) >= $k } $g->vertices;

    my @two_out = grep { $g->out_degree($_) == 2 &&
			 length $g->get_vertex_seq($_) >= $k } $g->vertices;

    my ($nin, $nout) = (scalar @two_in, scalar @two_out);
    warn "Combinatorial surprise: |two_in| = $nin, |two_out| = $nout"
	if $nin * $nout > 100;

    # Find undirected paths like  u<-v1->w  and  x->v3<-y  where one of
    # each of {u,w} and {x,y} is in common, and thus the two paths together
    # form the pattern we are looking for.

    my @quintets = ();
    foreach my $v1 (@two_out) {
	foreach my $v3 (@two_in) {
	    next if $v3 == $v1;

	    my ($u, $w, $x, $y) = ($g->successors($v1), $g->predecessors($v3));

	    # Swap each pair until if any two are in common, it is  w  and  x.
	    ($u, $w) = ($w, $u) if $u == $x || $u == $y;
	    ($x, $y) = ($y, $x) if $y == $w;

	    # Each instance of the pattern appears in the graph twice,
	    # so we pick the one with fewer negative vertices.
	    if ($w == $x && $u != $y && ! $g->has_edge($u, $y) &&
			count_negatives($u, $v1, $w, $v3, $y) <= 2) {
		push @quintets, [$u, $v1, $w, $v3, $y];
	    }
	}
    }

    return @quintets;
}

# Returns the vertices involved -- {u, w, x} -- if there is a crufty
# loop  v -> {u, w} -> x  hanging off v, or an empty list otherwise.
# The vertex given,  v,  should not have any other out-edges, so it
# should be selected so that its quintet-edges point towards it.
sub _cruft {
    my ($g, $v) = @_;

    return unless $g->out_degree($v) == 2;

    my ($u, $w) = $g->successors($v);
    return unless abs($g->get_vertex_contig_length($u) -
		      $g->get_vertex_contig_length($w)) <= 2;

    my @uu = $g->successors($u);
    my @ww = $g->successors($w);
    return unless scalar(@uu) == 1 && scalar(@ww) == 1 && $uu[0] == $ww[0];

    my $x = $uu[0];
    return unless $g->out_degree($x) == 0;

    return ($u, $w, $x);
}

=head2 annotate_bubbles

Blah

    $n = annotate_bubbles($graph, $max_indel);

=cut

sub _try_bubble {
    my ($g, $difference, $v, $w_v) = @_;

    return 0 unless $g->out_degree($v) == 2;

    my ($u, $w) = $g->successors($v);
    return 0 unless $g->in_degree($u) == 1 && $g->in_degree($w) == 1 &&
	abs($g->get_vertex_contig_length($u) -
	    $g->get_vertex_contig_length($w)) <= $difference;

    my @uu = $g->successors($u);
    my @ww = $g->successors($w);
    return 0 unless scalar(@uu) == 1 && scalar(@ww) == 1 && $uu[0] == $ww[0];

    my $x = $uu[0];
    return 0 unless $g->in_degree($x) == 2;

    # At this point,  v -> {u,w} -> x  is indeed a bubble.

    my $w_u = max($g->get_edge_weight($v, $u), $g->get_edge_weight($u, $x));
    my $w_w = max($g->get_edge_weight($v, $w), $g->get_edge_weight($w, $x));
    my $w_x = min($w_u, $w_w);
    $w_x = max($w_x, $w_v) if defined $w_v;

    $g->set_vertex_weight($x, $w_x);
    $g->set_vertex_weight(-$x, $w_x);
    return 1 + _try_bubble($g, $difference, $x, $w_x);
}

sub annotate_bubbles {
    my ($g, $difference) = @_;
    my $n = 0;
    $n += _try_bubble($g, $difference, $_, undef) foreach $g->source_vertices();
    return $n;
}

=head2 quintet_weights

Returns the maximum vertex weight and minimum edge weight within any quintet.
Blah.

    ($wb, $wq) = quintet_weights($graph, @quintets);

=cut

sub quintet_weights {
    my $g = shift;
    my $vw = 0;
    my $ew = 1000000;

    foreach (@_) {
	my @v = @$_;
	$vw = max($vw, map { $g->has_vertex_weight($_)?
			     $g->get_vertex_weight($_) : 0 } @v);
	$ew = min($ew, $g->get_edge_weight($v[1], $v[0]),
		       $g->get_edge_weight($v[1], $v[2]),
		       $g->get_edge_weight($v[2], $v[3]),
		       $g->get_edge_weight($v[4], $v[3]));
    }

    return ($vw, $ew);
}

=head2 quintet_score

    $score = quintet_score($graph);
    ($score, @quintets) = quintet_score($graph);

Returns a niceness score for the graph; in list context, also returns a list
of references to 5-element arrays, each containing an instance of the quintet
pattern.

A score of 100 indicates a perfect graph with five vertices forming a quintet.
Points are deducted for isolated vertices, cruft hanging off the quintet,
and major points for extra cycles and large-scale graph cruftiness.

=cut

sub quintet_score {
    my ($g) = @_;
    local $_;

    my (@isolated, @edged);
    foreach ($g->vertices) {
	next unless $_ > 0;
	if ($g->is_isolated_vertex($_)) { push @isolated, $_ }
	else { push @edged, $_ }
    }

    my @quintets = find_breakpoints($g);

    my $nisolated = scalar(@isolated);
    my $nquintets = scalar(@quintets);

    my $score;
    if ($nquintets == 0) {
	# If no quintets, score is capped number of isolated vertices.
	$score = $nisolated;
	$score = 15 if $score > 15;
    }
    elsif ($nquintets >= 2) {
	# Too many?  Score is middling, dropping with each extra quintet.
	# FIXME
	$score = 60 - 5 * ($nquintets - 2) - $nisolated;
    }
    else {
	$score = 100;
	my @q = @{$quintets[0]};
	print "Number of quintets: $nquintets\n" if $^D;
	print "Quintet: @q" if $^D;

	$g->set_vertex_attribute($_, 'style', 'filled') foreach @q;

	my %touched;
	$touched{abs $_}++ foreach @q;

	foreach my $v (0, -1, 3, -4) {
	    my @loop = _cruft($g, ($v >= 0)? $q[$v] : -$q[-$v]);
	    next unless @loop;

	    print ' with v', abs $v, " cruft @loop" if $^D;
	    $touched{abs $_}++ foreach @loop;
	    $g->set_vertex_attribute($_, 'style', 'filled') foreach @loop;
	    $g->set_vertex_attribute($loop[0], 'fillcolor', 'salmon');
	    $g->set_vertex_attribute($loop[1], 'fillcolor', 'salmon');
	    $g->set_vertex_attribute($loop[2], 'fillcolor', 'darksalmon');
	    $score -= 3;
	}
	print "\n" if $^D;

	# Big points off for any further crufty vertices
	foreach (@edged) {
	    $score -= 5 if !defined($touched{$_}) || $touched{$_} != 1;
	}

	# And points off for isolated vertices
	$score -= $nisolated;

	# And for smoothed cruft
	foreach (@q) {
	    my $cost = $g->get_vertex_attribute($_, 'bubble') || 0;
	    $score -= int(($cost + 1) / 2);
	}
    }
    print "Score: $score\n" if $^D;
    $score = 0 if $score < 0;
    return wantarray? ($score, @quintets) : $score;
}

=head2 concatenate_path

    $v = concatenate_path($graph, $u, $w, $v, [...]);

Merges the vertices along the path (and the twins on the twin path).
New vertex gets total of vertex weights and bubble attribute.
Returns the new vertex, or undef if it wasn't really a path.
OR SHOULD IT CROAK, LIKE IT DOES NOW?

=cut

sub concatenate_path {
    my $g = shift @_;

    my %attr;
    my $seq = '';
    foreach (@_) {
	$seq .= $g->get_vertex_seq($_);
	foreach my $key ('weight', 'cost', 'bubble') {
	    my $value = $g->get_vertex_attribute($_, $key);
	    $attr{$key} += $value if defined $value;
	}
    }
    my $twinseq = '';
    $twinseq .= $g->get_vertex_seq(-$_) foreach reverse @_;

    my $u = shift @_;

    my $p = $u;
    foreach (@_) {
	croak "$p -> $_ cannot be concatenated" unless $g->has_edge($p, $_)
			&& $g->out_degree($p) == 1 && $g->in_degree($_) == 1;
	$p = $_;
    }

    my $v = $_[-1];
    foreach ($g->successors($v)) {
	$g->add_weighted_edge($u, $_, $g->get_edge_weight($v, $_));
	$g->add_weighted_edge(-$_, -$u, $g->get_edge_weight(-$_, -$v));
    }

    foreach (@_) { $g->delete_vertices($_, -$_) }

    $g->set_vertex_attribute($u,  'seq', $seq);
    $g->set_vertex_attribute(-$u, 'seq', $twinseq);
    foreach my $key (keys %attr) {
	$g->set_vertex_attribute($u,  $key, $attr{$key});
	$g->set_vertex_attribute(-$u, $key, $attr{$key});
    }
    return $u;
}

=head2 try_smooth_bubble

    $u = try_smooth_bubble($graph, $u);

If  $u  is the base of a bubble (i.e., blah)
Pops the bubble  u -> {w1, w2} -> v  and its twin.  FIXME
Or returns undef if it wasn't a bubble.

=cut

sub try_smooth_bubble {
    my ($g, $u) = @_;

    my @w = $g->successors($u);
    return unless scalar(@w) >= 2;

    my @w0_successors = $g->successors($w[0]);
    return unless scalar(@w0_successors) == 1;
    my $v = $w0_successors[0];
    return unless $g->in_degree($v) == scalar(@w);

    foreach (@w) {
	return unless $g->in_degree($_) == 1 && $g->has_edge($_, $v)
		   && $g->out_degree($_) == 1;
    }

    # FIXME Cases involving twins are probably implementable, but punt for now.
    my %twins;
    foreach ($u, @w, $v) { return unless $twins{abs $_}++ == 0 }

    # At this point,  u -> {w0, w1, ...} -> v really is a bubble.

    my $w = update_bubble_consensus($g, @w);
    foreach (@w) { $g->delete_vertices($_, -$_) unless $_ == $w }
    return concatenate_path($g, $u, $w, $v);
}

=head2 update_bubble_consensus

    $w = update_bubble_consensus($graph, @w);

Modifies one of the bubble side-vertices (as listed in @w) and its twin, and
returns which pair was modified.  The vertices' C<seq> attributes and various
weights are modified to reflect the consensus of both/all sides of the bubble.

=cut

sub snp_cost { return 1; }
sub indel_cost { return $_[0]; }

sub update_bubble_consensus {
    my $g = shift @_;

    # TODO For now, the "consensus" is just the highest-coverage alternative.
    my $w = $_[0];
    my $w_cov = $g->get_vertex_weight($w);
    foreach (@_) {
	my $cov = $g->get_vertex_weight($_);
	if ($cov > $w_cov) { $w = $_; $w_cov = $cov; }
    }

    my $cost = 0;
    my $w_len = length $g->get_vertex_seq($w);
    foreach (@_) {
	next if $_ == $w;
	my $len_diff = abs ($w_len - length $g->get_vertex_seq($_));
	$cost += ($len_diff == 0)? snp_cost() : indel_cost($len_diff);
    }
    $g->set_vertex_attribute($w, 'bubble', $cost);

    return $w;
}

=head2 smooth_fringe_bubbles

    $count = smooth_fringe_bubbles($graph);

Blah.

=cut

sub smooth_fringe_bubbles {
    my ($g) = @_;
    my $n = 0;
    foreach my $v ($g->source_vertices()) {
	while (defined ($v = try_smooth_bubble($g, $v))) { $n++ }
    }
    return $n;
}

=head2 format_quintet

    $text = format_quintet(@quintet);

Returns a UTF-8 string containing the path through the graph representing
the feature of interest, shown as the vertex ids separated by Unicode arrows.

=cut

sub format_quintet {
    my $left  = "\x{2190}"; #  LEFTWARDS ARROW
    my $right = "\x{2192}"; # RIGHTWARDS ARROW

    return "$_[0] $left $_[1] $right $_[2] $right $_[3] $left $_[4]";
}

=head2 format_quintet_alignment

    $text = format_quintet_alignment($graph, @quintet);

Returns a 5-line string, formatted to show how the five contigs pile up
across the breakpoint.

=cut

sub format_quintet_alignment {
    my ($graph, @v) = @_;

    my $kminus1 = $graph->get_graph_hash_length() - 1;
    my @ctg = map $graph->get_vertex_contig($_), @v;
    my $lenctg2 = length($ctg[2]);

    my $overhang1 = length($ctg[1]) - $kminus1;
    my $overhang4 = length($ctg[4]) - $lenctg2;

    my ($overhangtab, $tab1, $tab4);
    if ($overhang1 > $overhang4) {
	$overhangtab = ' ' x $overhang1;
	$tab1 = '';
	$tab4 = ' ' x ($overhang1 - $overhang4);
    }
    else {
	$overhangtab = ' ' x $overhang4;
	$tab1 = ' ' x ($overhang4 - $overhang1);
	$tab4 = '';
    }

    my $extra = ' ' x ($lenctg2 - $kminus1);

    return <<"EOS"
$overhangtab$ctg[0]
$tab1$ctg[1]
$overhangtab$ctg[2]
$overhangtab$extra$ctg[3]
$tab4$ctg[4]
EOS
}

=head2 get_quintet_triplet

    ($one, $broken, $two) = get_quintet_triplet($graph, @quintet)

=cut

sub _concat {
    my ($graph, $v0, @v) = @_;
    my $contig = $graph->get_vertex_contig($v0);
    foreach (@v) { $contig .= $graph->get_vertex_seq($_) }
    return $contig;
}

sub get_quintet_triplet {
    my ($graph, @v) = @_;
    return (_concat($graph, $v[1], $v[0]), _concat($graph, @v[1..3]),
	    _concat($graph, $v[4], $v[3]));
}

=head2 get_bp_region

    ($within, $width) = get_bp_region($graph, @quintet);

Returns the sequence within the breakpoint.  If $within is empty and $width is
zero, it is a perfectly clean break.  If $width is positive, then $within is a
shard of non-templated sequence within the break.  S<If $width> is negative,
then $within is a micro-homology spanning the possible exact locations of the
break.

=cut

sub get_bp_region {
    my ($graph, @v) = @_;
    local $_;

    my $kminus1 = $graph->get_graph_hash_length() - 1;
    my $gap = $graph->get_vertex_contig_length($v[2]) - 2 * $kminus1;

    if ($gap >= 0) {
	$_ = substr($graph->get_vertex_seq($v[2]), 0, -$kminus1);
    }
    else {
	my $overlap = -$gap;
	$_ = substr($graph->get_vertex_contig($v[3]), 0, $overlap);
	return unless $_ eq substr($graph->get_vertex_contig($v[1]), -$overlap);
    }

    return wantarray? ($_, $gap) : $_;
}

=head2 get_bp_surrounding_region

    ($left, $within, $right, $width) =
	get_bp_surrounding_region($graph, @quintet);

Returns the sequence within the breakpoint, as per C<get_bp_region()>,
and also the sequences on either side of the breakpoint.

=cut

sub get_bp_surrounding_region {
    my ($graph, @v) = @_;

    my $ctg1 = $graph->get_vertex_contig($v[1]);
    my $ctg3 = $graph->get_vertex_contig($v[3]);

    my ($within, $gap) = get_bp_region($graph, @v);
    my $clip = ($gap >= 0)? 0 : -$gap;

    return (substr($ctg1, 0, length($ctg1) - $clip), $within,
	    substr($ctg3, $clip), $gap);
}

=head2 get_bp_alignment_text

    map_contigs($graph, @quintet); #FIXME
    $text = get_bp_alignment_text($graph, @quintet);

Returns the sequence within the breakpoint and the coordinates of the
adjoining ends of the contigs on either side.  The vertices in @quintet must
first have been annotated via C<map_contigs()>.

=over 4

=item C<Chr.1  123] ATGC [456  Chr.2>

Square brackets indicate a shard of non-templated sequence within the
breakpoint, or if it is a perfectly clean break, the two contigs will be
shown separated by C<][>.

=item C<Chr.1  119(23)--ATGC--456(60)  Chr.2>

Hyphens indicate that the two contigs overlap in a micro-homology, and that
the exact location of the breakpoint cannot be determined but is somewhere
in the range indicated.

=back

=cut

sub _get_alignment_summary {
    my ($graph, $v, $end) = @_;

    my $mappings = $graph->get_vertex_attribute($v, 'mappings');
    return ('???', undef) unless ref($mappings) eq 'ARRAY';
    return ('', undef) if scalar(@$mappings) == 0;
    my $aln = $$mappings[0];

    my $chr = $aln->target->prettyname;
    if ($end eq 'L') {
	my $qs = $aln->query->start;
	$chr .= "[\@$qs]" unless $qs == 1;
    }
    else {
	my $qe = $aln->query->end;
	$chr .= "[\@$qe]" unless $qe == $graph->get_vertex_contig_length($v);
    }

    my $n = scalar(@$mappings) - 1;
    if ($n > 0) {
	local $_ = ($n > 1)? 's' : '';
	$chr .= ' [score:' . $aln->score . ", $n other$_]";
    }

    return ($chr, $aln);
}

sub _compress_range {
    my $a = "$_[0]";
    my $b = "$_[1]";

    if (length($a) == length($b)) {
	return $a if $a eq $b;

	my $i = 0;
	$i++ while substr($a, $i, 1) eq substr($b, $i, 1);

	# Write a single digit difference as 216(17) rather than 216(7).
	$i-- if $i > 0 && $i == length($b) - 1;

	$b = substr($b, $i);
    }

    return "$a($b)";
}

sub _breakpoint_summary {
    my ($text, $gap, $left, $right) = @_;

    my $clip = ($gap >= 0)? 0 : -$gap;
    my ($L, $R);

    if (defined $left) {
	my $clipleft = $left->target->clone()->trim(0, $clip);
	$L = _compress_range($clipleft->end, $left->target->end);
#$L .='{'._compress_range($left->target->end - $clip, $left->target->end).'}';
    }
    else { $L = 'unmappable' }

    if (defined $right) {
	my $clipright = $right->target->clone()->trim($clip, 0);
	$R = _compress_range($right->target->start, $clipright->start);
#$R .='{'._compress_range($right->target->start, $right->target->start + $clip).'}';
    }
    else { $R = 'unmappable' }

    substr($text, 17) = "..." if length($text) > 20; # FIXME NUKE-ME

    if ($gap > 0)    { return "  $L] $text [$R  " }
    elsif ($gap < 0) { return "  $L--$text--$R  " }
    else             { return "  $L][$R  " }
}

sub get_bp_alignment_text {
    my ($graph, @v) = @_;

    my ($name1, $aln1) = _get_alignment_summary($graph, $v[1], 'R');
    my ($name3, $aln3) = _get_alignment_summary($graph, $v[3], 'L');

    my ($text, $gap) = get_bp_region($graph, @v);

    my $category = (defined $aln1 && defined $aln3 &&
		    $aln1->target->name   eq $aln3->target->name &&
		    $aln1->target->strand eq $aln3->target->strand)?
		   'indel' : 'rearrangement';

    local $_ = $name1 . _breakpoint_summary($text, $gap, $aln1, $aln3) . $name3;
    return wantarray? ($_, $category) : $_;
}

=head2 get_bp_alignment

    ($chr1, $pretty1, $pos1, $pos1h, $strand1, $seq, $width,
     $chr2, $pretty2, $pos2, $pos2h, $strand2) = get_bp_alignment($graph, @quintet);

=cut

sub get_bp_alignment {
    my ($graph, @v) = @_;

    my ($name1, $aln1) = _get_alignment_summary($graph, $v[1], 'R');
    my ($name3, $aln3) = _get_alignment_summary($graph, $v[3], 'L');

    my ($text, $gap) = get_bp_region($graph, @v);
    return unless defined $gap;
    my $clip = ($gap >= 0)? 0 : -$gap;

    my @ret;
    if (defined $aln1) {
	my $clipleft = $aln1->target->clone()->trim(0, $clip);
	push @ret, $aln1->target->name, $name1, $clipleft->end, $aln1->target->end, $aln1->target->strand;
    }
    else { push @ret, undef, undef, 0, 0, undef; }

    push @ret, $text, $gap;

    if (defined $aln3) {
	my $clipright = $aln3->target->clone()->trim($clip, 0);
	push @ret, $aln3->target->name, $name3, $aln3->target->start, $clipright->start, $aln3->target->strand;
    }
    else { push @ret, undef, undef, 0, 0, undef; }

    return @ret;
}

=head2 get_isolated_bp_alignment

    $description = get_isolated_bp_alignment(@mappings);

Returns the sequence within the breakpoint and the coordinates of the
adjoining ends of the contigs on either side.  Of a breakpoint in both
copies, as represented by the (sorted) mappings of a single vertex.
Or C<undef>, if the mappings do not constitute a set of breakpoints.

=cut

sub _breakpoint_summary_iso {
    my ($text, $gap, $left, $right) = @_;
    my ($tleft, $tright) = ($left->target, $right->target);

    substr($text, 17) = '...' if length($text) > 20; # FIXME NUKE-ME

    if ($gap >= 0) {
	$text = " $text " if $gap > 0;
	return '  '. $tleft->end ."]$text\[". $tright->start .'  ';
    }
    else {
	my $overlap = -$gap;
	my $tleft2  =  $tleft->clone()->trim(0, $overlap);
	my $tright2 = $tright->clone()->trim($overlap, 0);
	my $L = _compress_range($tleft2->end, $tleft->end);
	my $R = _compress_range($tright->start, $tright2->start);

#	$L .='{'._compress_range($tleft->end - $overlap, $tleft->end).'}';
#	$R .='{'._compress_range($tright->start, $tright->start + $overlap).'}';
	return "  $L--$text--$R  ";
    }
}

sub get_isolated_bp_alignment {
    my ($graph, $v) = @_;

    my $mappings = $graph->get_vertex_attribute($v, 'mappings');
    return undef if ref($mappings) ne 'ARRAY';
    my @mappings = @$mappings;
    return undef if scalar(@mappings) < 2;

    my $prevprevend = 0;
    my $prevend = 0;
    foreach (@mappings) {
	my $q = $_->query;
	return undef if $q->end <= $prevend || $q->start <= $prevprevend;
	$prevprevend = $prevend;
	$prevend = $q->end;
    }

    my $query_text = $graph->get_vertex_contig($v);

    # The output is returned in the form CHRNAME[  END SHARD START  CHRNAME]*
    my $output = "";
    my %targets;

    my $left = undef;
    foreach my $right (@mappings) {
	if (defined $left) {
	    my $rqstart = $right->query->start;
	    my $lqlim = $left->query->end + 1;
	    my $gap = $rqstart - $lqlim;
	    my $text = substr($query_text, min($lqlim,$rqstart) - 1, abs($gap));
	    $output .= _breakpoint_summary_iso($text, $gap, $left, $right);
	}

	$output .= $right->target->prettyname;
	$targets{$right->target->name.$right->target->strand}++;
	$left = $right;
    }

    my $category = (scalar(keys %targets) > 1)? 'rearrangement' : 'indel';
    return wantarray? ($output, $category) : $output;
}

sub get_isolated_bp_surrounding_region {
    my ($graph, $v) = @_;

    my $mappings = $graph->get_vertex_attribute($v, 'mappings');
    return '???' unless ref($mappings) eq 'ARRAY';
    my @mappings = @$mappings;
    return 'unmappable' if scalar(@mappings) == 0;
    return undef if scalar(@mappings) < 2;

    my $query_text = $graph->get_vertex_contig($v);
    my $left = $mappings[0];
    my $right = $mappings[1];

    my $rqstart = $right->query->start;
    my $lqlim = $left->query->end + 1;
    my $gap = $rqstart - $lqlim;
    my $within = substr($query_text, min($lqlim,$rqstart) - 1, abs($gap));

    my $clip = ($gap >= 0)? 0 : -$gap;
    $left->query->trim(0, $clip);
    $right->query->trim($clip, 0);

    return ($left->query->seq, $within, $right->query->seq, $gap);
}

1;
__END__

=head2 trim_alignment

    $subsugar = trim_alignment($seqtext, $subsugar, $refseqtext);
    $subsugar = trim_alignment($st, $ss, $rst, $ends, $threshold);

Blah blah.

=cut

sub _trim_end {
    my ($query, $qbase, $target, $tbase, $dir, $len, $threshold) = @_;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
			     eq substr($target, $tbase + $dir * $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query,  $qbase + $dir * $i, 1)
			     ne substr($target, $tbase + $dir * $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $threshold;

    return $startpos;
}

sub trim_alignment {
    my ($query, $subsugar, $subtarget, $ends, $threshold) = @_;
    $ends = 'LR' unless defined $ends;
    $threshold = 6 unless defined $threshold;

    my ($qpos, $qlim, $qstrand, $tname, $tpos, $tlim, $tstrand, $score) =
	split /\s+/, $subsugar;

    if ($ends =~ /[Ll]/) {
	my $trim = _trim_end($query, $qpos, $subtarget, 0, +1,
			     length $subtarget, $threshold);
	$qpos += $trim;
	$tpos += $trim;
    }

    if ($ends =~ /[Rr]/) {
	my $trim = _trim_end($query, $qlim - 1, $subtarget, -1, -1,
			     length $subtarget, $threshold);
	$qlim -= $trim;
	$tlim -= $trim;
    }

    return "$qpos $qlim $qstrand $tname $tpos $tlim $tstrand $score";
}

=head2 trim_alignment

    Bio::Brass::trim_threshold = 6;
    @subsugar = trim_alignment_left ($seqtext, $refseqtext, @subsugar);
    @subsugar = trim_alignment_right($seqtext, $refseqtext, @subsugar);
    @subsugar = trim_alignment      ($seqtext, $refseqtext, @subsugar);

Blah blah.

=cut

our $trim_threshold = 6;

use constant { QPOS => 0, QLIM => 1, TPOS => 4, TLIM => 5 };

# Subsugar is  qpos qlim qstrand tname tpos tlim tstrand score


# FIXME  In all this, adjusting SUBSUGAR means that SUBTARGET is now too big;
# we possibly ought to trim it too.

sub trim_alignment_left {
    my ($query, $subtarget, @subsugar) = @_;
    my $qpos = $subsugar[QPOS];
    my $len = length $subtarget;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query, $qpos + $i, 1)
			     eq substr($subtarget, $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query, $qpos + $i, 1)
			     ne substr($subtarget, $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $trim_threshold;

    $subsugar[QPOS] += startpos;
    $subsugar[TPOS] += startpos;
    return @subsugar;
}

#	my $trim = _trim_end($query, $qlim - 1, $subtarget, -1, -1,
#			     length $subtarget, $threshold);
sub trim_alignment_right {
    my ($query, $subtarget, @subsugar) = @_;
    my $qend = $subsugar[QLIM] - 1;
    my $len = length $subtarget;
    my ($startpos, $overhang);

    my $i = 0;
    do {
	$startpos = $i;

	$i++ while $i < $len && substr($query, $qend - $i, 1)
			     eq substr($subtarget, -1 - $i, 1);
	my $matched_count = $i - $startpos;
	my $mismatchpos = $i;
	$i++ while $i < $len && substr($query, $qend - $i, 1)
			     ne substr($subtarget, -1 - $i, 1);
	my $mismatched_count = $i - $mismatchpos;
	$overhang = $matched_count - $mismatched_count;
    } until $i == $len || $overhang >= $threshold;

    $subsugar[QLIM] -= startpos;
    $subsugar[TLIM] -= startpos;
    return @subsugar;
}

sub trim_alignment {
    my ($query, $subtarget, @subsugar) = @_;

    @subsugar = trim_alignment_left ($query, $subtarget, @subsugar);
    @subsugar = trim_alignment_right($query, $subtarget, @subsugar);
    return @subsugar;

    return trim_alignment_left($query, $subtarget,
		    trim_alignment_right($query, $subtarget, @subsugar));
}

=head1 SEE ALSO

=over 4

=item L<Graph>

Jarkko Hietaniemi's classes for representing graphs.

=back

=head1 AUTHOR

John Marshall E<lt>jm18@sanger.ac.ukE<gt>

=cut

1;
