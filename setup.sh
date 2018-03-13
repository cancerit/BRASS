#!/bin/bash

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


SOURCE_BLAT="http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"

# if issues found downgrade to 2.23.0 but can't find any use of bedtools coverage
SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz"

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

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path  such as /opt/pancan"
  exit 0
fi

INST_PATH=$1

if [[ "x$2" == "x" ]] ; then
  INST_METHOD=0
else
  INST_METHOD=$2
fi

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

set -eu

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
PERLROOT=$INST_PATH/lib/perl5

# allows user to knowingly specify other PERL5LIB areas.
if [ -z ${CGP_PERLLIBS+x} ]; then
  export PERL5LIB="$PERLROOT"
else
  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
fi

#add bin path for install tests
export PATH=$INST_PATH/bin:$PATH

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`
echo $CPANM

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

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Vcf`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install cgpVcf before proceeding: https://github.com/cancerit/cgpVcf/releases"
  exit 1;
fi

perlmods=( "Graph" )

set +e
for i in "${perlmods[@]}" ; do
  echo "Installing build prerequisite $i..."
  $CPANM --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

cd $SETUP_DIR

set -e

if [ $INST_METHOD -eq 1 ] ; then
  echo -e "\n\t !!! Not installing additional tools INST_METHOD requested !!! \n\n"
else
  echo "Building bedtools2 ..."
  if [ -e $SETUP_DIR/bedtools.success ]; then
    echo " previously installed ...";
  else
    cd $SETUP_DIR
    get_distro "bedtools2" $SOURCE_BEDTOOLS
    mkdir -p bedtools2
    tar --strip-components 1 -C bedtools2 -zxf bedtools2.tar.gz
    make -C bedtools2 -j$CPU
    cp bedtools2/bin/* $INST_PATH/bin/.
    touch $SETUP_DIR/bedtools.success
  fi


  echo "Building blat ..."
  if [ -e $SETUP_DIR/blat.success ]; then
    echo " previously installed ..."
  else
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
  fi


  cd $INIT_DIR
  echo "Building velvet..."
  if [ -e $SETUP_DIR/velvet.success ]; then
    echo " previously installed ..."
  else
    cd $INIT_DIR/distros
    tar -mzxf velvet_1.2.10.tgz
    cd velvet_1.2.10
    make MAXKMERLENGTH=95 velveth velvetg
    mv velveth $INST_PATH/bin/velvet95h
    mv velvetg $INST_PATH/bin/velvet95g
    make clean
    make velveth velvetg   	# don't do multi-threaded make
    mv velveth $INST_PATH/bin/velvet31h
    mv velvetg $INST_PATH/bin/velvet31g
    ln -fs $INST_PATH/bin/velvet95h $INST_PATH/bin/velveth
    ln -fs $INST_PATH/bin/velvet95g $INST_PATH/bin/velvetg
    cd $INIT_DIR
    rm -rf $INIT_DIR/distros/velvet_1.2.10
    touch $SETUP_DIR/velvet.success
  fi

  cd $INIT_DIR
  echo "Building exonerate..."
  if [ -e $SETUP_DIR/exonerate.success ]; then
    echo " previously installed ..."
  elif [ $INST_METHOD == 2 ]; then
    echo " Skipping exonerate install ..."
  else
    cd $INIT_DIR/distros
    tar zxf exonerate-2.2.0.tar.gz
    cd $INIT_DIR/distros/exonerate-2.2.0
    cp $INIT_DIR/distros/patches/exonerate_pthread-asneeded.diff .
    patch -p1 < exonerate_pthread-asneeded.diff
    ./configure --prefix=$INST_PATH
    make    # don't do multi-threaded make
    make check
    make install
    cd $INIT_DIR
    rm -rf $INIT_DIR/distros/exonerate-2.2.0
    touch $SETUP_DIR/exonerate.success
  fi
fi

cd $INIT_DIR
echo "Building brass (c++)..."
if [ -e $SETUP_DIR/brass.success ]; then
  echo " previously installed ..."
else
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
fi

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR/perl

echo "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
$CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null

echo "Installing brass (perl)..."
cd $INIT_DIR/perl
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
