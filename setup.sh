#!/bin/bash

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


SOURCE_BLAT="http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"
SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.21.0/bedtools-2.21.0.tar.gz"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/pancan"
  exit 0
fi

INST_PATH=$1

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

# get current directory
INIT_DIR=`pwd`

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

set -eu

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
PERLARCH=$PERLROOT/$ARCHNAME
export PERL5LIB="$PERLROOT:$PERLARCH"

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' PCAP`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install PCAP-core before proceeding: https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
  exit 1;
fi

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Grass`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install grass before proceeding: https://github.com/cancerit/grass/releases"
  exit 1;
fi

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

## grab cpanm:
rm -f $SETUP_DIR/cpanm
get_file $SETUP_DIR/cpanm http://xrl.us/cpanm
chmod +x $SETUP_DIR/cpanm

perlmods=( "Graph" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    perl $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo
  ) >>$INIT_DIR/setup.log 2>&1
  done_message "" "Failed during installation of $i."
done

cd $SETUP_DIR

echo -n "Building bedtools ..."
if [ -e $SETUP_DIR/bedtools.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  (
  set -x
  if [ ! -e bedtools ]; then
    get_distro "bedtools2" $SOURCE_BEDTOOLS
    mkdir -p bedtools2
    tar --strip-components 1 -C bedtools2 -zxf bedtools2.tar.gz
  fi
  make -C bedtools2 -j$CPU
  cp bedtools2/bin/* $INST_PATH/bin/.
  touch $SETUP_DIR/bedtools.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bedtools."

echo -n "Building blat ..."
if [ -e $SETUP_DIR/blat.success ]; then
  echo -n " previously installed ..."
else
  (
    get_distro "blat" $SOURCE_BLAT
    unzip -qu blat.zip
    cd $SETUP_DIR/blatSrc
    BINDIR=$SETUP_DIR/blat/bin
    export BINDIR
    export MACHTYPE
    mkdir -p $BINDIR
    make -j$CPU
    cp $BINDIR/blat $INST_PATH/bin/.
    touch $SETUP_DIR/blat.success
  ) >>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build blat."

cd $INIT_DIR
echo -n "Building brass (c++)..."
if [ -e $SETUP_DIR/brass.success ]; then
  echo -n " previously installed ..."
else
  (
    rm -rf $INIT_DIR/cansam*
    unzip -q distros/cansam.zip
    mv cansam-master cansam
    make -C cansam
    make -C c++
    cp c++/augment-bam $INST_PATH/bin/.
    cp c++/brass-group $INST_PATH/bin/.
    cp c++/filterout-bam $INST_PATH/bin/.
    make -C c++ clean
    rm -rf cansam
    touch $SETUP_DIR/brass.success
  ) >>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build brass (c++)."

cd $INIT_DIR
echo -n "Building velvet..."
if [ -e $SETUP_DIR/velvet.success ]; then
  echo -n " previously installed ..."
else
  (
    cd $INIT_DIR/distros
    tar zxf velvet_1.2.10.tgz
    cd velvet_1.2.10
    make MAXKMERLENGTH=95 velveth velvetg
    mv velveth $INST_PATH/bin/velvet95h
  	mv velvetg $INST_PATH/bin/velvet95g
  	make clean
  	# don't do multi-threaded make
  	make velveth velvetg

  	mv velveth $INST_PATH/bin/velvet31h
  	mv velvetg $INST_PATH/bin/velvet31g
  	ln -fs $INST_PATH/bin/velvet95h $INST_PATH/bin/velveth
  	ln -fs $INST_PATH/bin/velvet95g $INST_PATH/bin/velvetg

  	cd $INIT_DIR
  	rm -rf $INIT_DIR/distros/velvet_1.2.10

    touch $SETUP_DIR/velvet.success
  ) >>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build velvet."

cd $INIT_DIR
echo -n "Building exonerate..."
if [ -e $SETUP_DIR/exonerate.success ]; then
  echo -n " previously installed ..."
else
  (
    cd $INIT_DIR/distros
    tar zxf exonerate-2.2.0.tar.gz
    cd $INIT_DIR/distros/exonerate-2.2.0
    cp $INIT_DIR/distros/patches/exonerate_pthread-asneeded.diff .
    patch -p1 < exonerate_pthread-asneeded.diff
    ./configure --prefix=$INST_PATH
    # don't do multi-threaded make
    make
    make check
    make install;

    cd $INIT_DIR
    rm -rf $INIT_DIR/distros/exonerate-2.2.0

    touch $SETUP_DIR/exonerate.success
  ) >>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build exonerate."

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR/perl

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  perl $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

echo -n "Installing brass (perl)..."
(
  cd $INIT_DIR/perl
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
) >>$INIT_DIR/setup.log 2>&1
done_message "" "brass (perl) install failed."

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
