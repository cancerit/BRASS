# Makefile for Brass package, which identifies breakpoints via local assembly.

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

CANSAM = ../cansam

CXXFLAGS = -Wall -Wextra -g -I$(CANSAM) -O2
LDFLAGS  = -L$(CANSAM)

ALL_TARGETS = augment-bam brass-group filterout-bam
all: $(ALL_TARGETS)

ALL_OBJS = $(AUGMENT_BAM_OBJS) $(BRASS_GROUP_OBJS) $(FILTEROUT_BAM_OBJS)

AUGMENT_BAM_OBJS = augment-bam.o
augment-bam: $(AUGMENT_BAM_OBJS) $(CANSAM)/libcansam.a
	$(CXX) $(LDFLAGS) -o $@ $(AUGMENT_BAM_OBJS) -lcansam -lz

BRASS_GROUP_OBJS = brass-group.o feature.o rearrgroup.o
brass-group: $(BRASS_GROUP_OBJS) $(CANSAM)/libcansam.a
#	$(CXX) $(LDFLAGS) -o $@ $(BRASS_GROUP_OBJS) -lcansam -lz
	$(CXX) $(LDFLAGS) -o $@ $(BRASS_GROUP_OBJS) -lcansam -lboost_iostreams -lz

FILTEROUT_BAM_OBJS = filterout-bam.o
filterout-bam: $(FILTEROUT_BAM_OBJS) $(CANSAM)/libcansam.a
	$(CXX) $(LDFLAGS) -o $@ $(FILTEROUT_BAM_OBJS) -lcansam -lz

ALL_TEMPORARIES = version.inc TAGS

augment-bam.o: augment-bam.cpp version.inc
brass-group.o: brass-group.cpp feature.h imergestream.h rearrgroup.h \
	       version.inc
feature.o: feature.cpp feature.h
filterout-bam.o: filterout-bam.cpp version.inc
rearrgroup.o: rearrgroup.cpp rearrgroup.h

version.inc: ../perl/lib/Bio/Brass.pm
	perl -I../perl/lib -MBio::Brass -e 'print <<"EOF"'		\
	  -e '/* Generated from lib/Bio/Brass.pm */'		\
	  -e '#define BRASS_VERSION "$VERSION"'	\
	  -e 'EOF' > $@


.PHONY: all clean tags install

clean:
	-rm -f $(ALL_TARGETS) $(ALL_OBJS) $(ALL_TEMPORARIES)

tags:
	ctags -o TAGS *.h *.cpp lib/*/*.pm lib/*/*/*.pm

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
datarootdir = $(prefix)/share
mandir      = $(datarootdir)/man
man1dir     = $(mandir)/man1

INSTALL = install
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644
install: $(ALL_TARGETS)
	mkdir -p $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
	rsync -av --chmod=Du=rwx,Dg=rx,Do=rx,Fu=rw,Fg=r,Fo=r --exclude=".*" lib $(DESTDIR)$(bindir)/
	$(INSTALL_PROGRAM) augment-bam $(DESTDIR)$(bindir)/augment-bam
	$(INSTALL_PROGRAM) brass-group $(DESTDIR)$(bindir)/brass-group
	$(INSTALL_PROGRAM) filterout-bam $(DESTDIR)$(bindir)/filterout-bam
	$(INSTALL_PROGRAM) group-rearrs $(DESTDIR)$(bindir)/group-rearrs
	$(INSTALL_PROGRAM) make-repeat-file $(DESTDIR)$(bindir)/make-repeat-file
	$(INSTALL_PROGRAM) rearr2bedpe $(DESTDIR)$(bindir)/rearr2bedpe
	$(INSTALL_PROGRAM) brass-assemble $(DESTDIR)$(bindir)/brass-assemble
	$(INSTALL_DATA) brass-group.1 $(DESTDIR)$(man1dir)/brass-group.1
	$(INSTALL_DATA) rearr2bedpe.1 $(DESTDIR)$(man1dir)/rearr2bedpe.1
	$(INSTALL_DATA) brass-assemble.1 $(DESTDIR)$(man1dir)/brass-assemble.1
