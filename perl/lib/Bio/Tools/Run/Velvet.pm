package Bio::Tools::Run::Velvet;

########## LICENCE ##########
# Copyright (c) 2014,2015 Genome Research Ltd.
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

use Bio::Root::IO;
use Bio::SeqIO;

use base qw(Bio::Root::Root
	    Bio::Tools::Run::WrapperBase
	    Bio::Factory::ApplicationFactoryI);

our $VERSION = '0.60'; # FIXME Remove for BioPerl submission

=head1 NAME

Bio::Tools::Run::Velvet - de novo assembly with Velvet

=head1 SYNOPSIS

    use Bio::Tools::Run::Velvet;

    $factory = Bio::Tools::Run::Velvet->new(@parameters);

    $factory->prepare(@sequence_files);

    $factory->run(@arguments);
    $factory->run(@more_arguments);
    ...

=head1 DESCRIPTION

B<Bio::Tools::Run::Velvet> is a wrapper module for I<Velvet>,
a I<de novo> genomic assembler designed for short read sequencing technologies.

This module invokes the B<velvetg> and B<velveth> programs, which ideally will
be found on your path with no further ado.  See below for how to check this and
what to do if they are not installed on your path.

=head2 new

    $factory = Bio::Tools::Run::Velvet->new(-hash_length => 31);

Creates and returns a new Velvet wrapper object.  Any number of the following
named arguments may be given, and C<-hash_length> must be given.

=over 4

=item -hash_length => I<integer>

The hash length (or I<k>-mer length) to be used, which must be an
S<odd integer>.  This argument must be given; there is no default.

=item -dir => I<directory>

blah blah.

=item -verbose => I<integer>

blah blah.

=item -program_dir => I<directory>

=item -program_name => I<pattern>

Program directory and name pattern, equivalent to setting them via the
methods described below.

=back

=cut

sub new {
    my $class = shift @_;
    my $self = $class->SUPER::new(@_);

    $self->_set_from_args(\@_, -methods =>
	    [qw(program_dir program_name quiet save_tempfiles verbose)]);

    my ($k, $dir) = $self->_rearrange([qw(HASH_LENGTH DIR)], @_);

    if (defined $k) { $self->{'hash_length'} = $k }
    else { $self->throw('No -hash_length specified') }

    if (defined $dir) {
	# FIXME do we want to save other tempfiles too?
	# Or do we really want our own _save_working_dir property?
	$self->save_tempfiles(1);
	$self->tempdir($dir);
    }

    return $self;
}

=head2 prepare

    $factory->prepare('reads.fa', $seq, '-shortPaired', 'reads.fastq',
                      -eland => 'foobar', ...);

Runs B<velveth> to prepare the working directory, building the
Velvet data files corresponding to the sequence files specified.
Any number of sequence file arguments may be given, possibly interspersed
with read category options.

Each item in the argument list may be either:

=over 2

=item *

a file name; the file's format is determined by examining the filename
extension

=item *

two list items, -I<format> =E<gt> F<filename>, specifying a file and its
S<format (as> a B<velveth> file format option), for cases when this cannot
be determined from the extension

=item *

a Bio::Seq object, or a reference to an array of Bio::Seq objects

=item *

a B<velveth> read category option (C<-shortPaired>, C<-long>, etc), which
applies to subsequent files and sequences in the argument list.

=back

Returns the total number of sequences read from all the files specified,
or throws an exception if B<velveth> invocation failed.

=cut

sub prepare {
    my $self = shift @_;
    local $_;

    my $velveth = $self->executable('velveth');
    my $cmd = $velveth . q{ } . $self->tempdir() . q{ } . $self->{'hash_length'} . q{ -create_binary};

    my $short_counter = '';
    while (scalar(@_) > 0) {
	if (ref $_[0]) {
	    my $val = shift @_;

	    my ($TMP, $tmpfname) = $self->io->tempfile();
	    my $tmpio = Bio::SeqIO->new(-fh => $TMP, -format => 'fasta');

	    # Accept either an array reference or a single Bio::Seq object; if
	    # it's some other random reference, we rely on fasta's write_seq()
	    # to throw a suitable exception.
	    if (ref($val) eq 'ARRAY') { $tmpio->write_seq($_) foreach @$val }
	    else { $tmpio->write_seq($val) }

	    $cmd .= " -fasta -short$short_counter $tmpfname";
	}
	elsif ($_[0] =~ /^-(short|long)/) {
	    my $type = shift @_;
	    $cmd .= " $type";
	}
	elsif ($_[0] =~ /^-/) {
	    my $format = shift @_;
	    my $fname = shift @_;
	    $cmd .= " $format $fname";
	}
	else {
	    my $fname = shift @_;
	    $_ = lc $fname;
	    if (/\.f(ast)?([aq])((.gz)?)$/) { $cmd .= " -fast$2$3 -short$short_counter $fname" }
	    elsif (/\.sam((.gz)?)$/) { $cmd .= " -sam$1 -short$short_counter $fname" }
	    elsif (/\.bam$/) { $cmd .= " -bam -short$short_counter $fname" }
	    else { $self->throw("Can't determine format for $fname") }
	}
	$short_counter = 1 unless $short_counter;
	$short_counter++;
    }

    my $count = undef;

    # FIXME Output buffering

    print "$cmd\n" if $self->verbose >= 0;
    open my $VELVET, "$cmd|" or $self->throw("Can't execute $velveth: $!");
    while (<$VELVET>) {
	print unless $self->quiet || $self->verbose < 0;
	$count = $1 if /^(\d+)\s+sequences in total/;
    }
    close $VELVET or $self->throw("velveth execution failed");

    return $count;
}

=head2 run

    $factory->run(-cov_cutoff => 1.1, -read_trkg => 'yes', ...);

Runs B<velvetg> blah blah blah.

=cut

sub run {
    my $self = shift @_;
    local $_;

    my $velvetg = $self->executable('velvetg');
    my $cmd = join(" ", $velvetg,$self->tempdir(), '-shortMatePaired yes', '-read_trkg yes', '-unused_reads yes', @_);

    my $warnings = "";
    my $finalline;

    print "$cmd\n" if $self->verbose >= 0;
    open my $VELVET, "$cmd|" or $self->throw("Can't execute $velvetg: $!");
    while (<$VELVET>) {
	print unless $self->quiet || $self->verbose < 0;
    	$warnings .= $_ if substr($_, 0, 8) eq "WARNING:";
	$finalline = $_ if substr($_, 0, 12) eq "Final graph ";
    }
    close $VELVET or $self->throw("velvetg execution failed");

    $self->error_string($warnings);
    if (defined $finalline && $finalline =~
	    /(\d+)\s*nodes.*n50\s+of\s*(\d+).*max\s*(\d+).*total\s*(\d+)/) {
	return wantarray? ($1, $2, $3, $4) : $1;
    }

    return;
}

=head2 program_dir

    $dir = $factory->program_dir();

Returns the directory containing the B<velvetg> and B<velveth> executables,
as set by a previous invocation with an argument or by C<VELVETDIR> or
C<BIOTOOLDIR> environment variables, or returns C<undef>, indicating that
the executables are to be found by a path search.

    $factory->program_dir('/usr/local/bin');

With an argument, sets a program directory (and returns it), overriding the
normal setting that would otherwise apply for C<$factory>.

=cut

sub program_dir {
    my ($self, $dir) = @_;
    $self->{velvet_dir} = $dir if defined $dir;

    if (exists $self->{velvet_dir}) { return $self->{velvet_dir} }
    elsif (exists $ENV{VELVETDIR})  { return $ENV{VELVETDIR} }
    elsif (exists $ENV{BIOTOOLDIR}) { return $ENV{BIOTOOLDIR} }
    else { return undef }
}

=head2 program_name

    $name = $factory->program_name();

Returns a pattern matching the names of the Velvet executables, normally
C<velvet[gh]>.  This is used by C<executable()> to construct the names of
the particular B<velvetg> and B<velveth> executables to be invoked.

    $factory->program_name('velvet[gh]_de');

With an argument, sets the name pattern (and returns it).  The example shown
tells C<$factory> to use the colorspace version of Velvet.

=cut

sub program_name {
    my ($self, $name) = @_;
    $self->{velvet_name} = $name if defined $name;
    return $self->{velvet_name} || 'velvet[gh]';
}

=head2 executable

    $filename = $factory->executable('velvetg');

Returns blah blah blah.

=cut

sub executable {
    my ($self, $name, $warn) = @_;
    $self->throw('No name specified') unless defined $name;

    local $_;
    if ($name =~ /^velvet([a-z])$/) {
	my $gh = $1;
	$_ = $self->program_name();
	s/\[[a-z]+\]/$gh/; # Replace [gh] and similar
    }
    else {
	$_ = $name;
    }

    my $dir = $self->program_dir();
    return ($dir)? Bio::Root::IO->catfile($dir, $_) : $_;
}

=head2 version

    $version = $factory->version();

Returns a string containing the version number of the B<velvetg> and B<velveth>
programs (or S<"I<v1>[velvetg] I<v2>[velveth]"> if the versions differ),
or C<undef> if an executable cannot be found or a version number is
otherwise unavailable.

=cut

sub _run_for_version {
    my ($self, $executable) = @_;
    local $_;

    my $FH;
    if ($self->io->exists_exe($executable) && open $FH, "-|", $executable) {
	while (<$FH>) { return $1 if /^Version\s+(\S+)/ }
    }

    return undef;
}

sub version {
    my ($self) = @_;

    my $velvetg = $self->executable('velvetg');
    my $velveth = $self->executable('velveth');
    my $g = $self->_run_for_version($velvetg);
    my $h = $self->_run_for_version($velveth);
    return undef unless defined $g && defined $h;

    return ($g eq $h)? $g : "${g}[$velvetg] ${h}[$velveth]";
}

=head2 outfile

Blah.

=cut

sub outfile {
    my ($self, $basename) = @_;
    return Bio::Root::IO->catfile($self->tempdir(), $basename);
}

1;
__END__

=head2 Helping the module find your executables

If typing B<which velvetg> and B<which velveth> on your system displays the
full pathnames of suitable executables, then this module will also be able
to find these executables via this path.

Otherwise, or if you want to override the executables found via the path,
you can set an environment variable C<VELVETDIR> to the full pathname of
a directory containing the B<velvetg> and B<velveth> executables.
(Or, if that variable is unset, C<BIOTOOLDIR> is consulted similarly.)

In sh or bash
    export VELVETDIR=/home/username/build/velvet_x.y.z

In csh or tcsh
    setenv VELVETDIR /home/username/build/velvet_x.y.z

Or hard-code the directory into each script that uses this module:

    use Bio::Tools::Run::Velvet;
    $ENV{VELVETDIR} = '/home/username/build/velvet_x.y.z';

If you are running a script elsewhere, such as on a webserver, ensure that
B<that> environment also has a suitable path or sets a suitable environment
variable using one of these alternatives.

In detail, the executables are found: in the directory named by a C<VELVETDIR>
or (if that is unset) a C<BIOTOOLDIR> environment variable; or, if neither
of those is set, via the usual path search.

=head1 SEE ALSO

=over 4

=item L<http://www.ebi.ac.uk/~zerbino/velvet/>

The Velvet home page, with links to source code, documentation, papers,
and mailing list.

=item L<Graph::Reader::Velvet>

A class for reading the F<LastGraph> file produced by Velvet into an instance
of Jarkko Hietaniemi's L<Graph> class.

=back

=head1 AUTHOR

John Marshall E<lt>john.marshall@sanger.ac.ukE<gt>

=cut
