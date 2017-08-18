#!/bin/bash

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

if [ "$#" -lt "1" ] ; then
  echo "USAGE: ./setupR.sh /install/path"
  echo -e "\tOptionally add '1' to command to request build of R source in the install location"
  exit 1
fi

INST_PATH=$1
BUILD_R=$2

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

set -ue

# get current directory
INIT_DIR=`pwd`

#add bin path for use in R lib building
export PATH=$INST_PATH/bin:$PATH

# set R lib paths based on INST_PATH
export R_LIBS=$INST_PATH/R-lib
export R_LIBS_USER=$INST_PATH/R-lib

mkdir -p $R_LIBS

TMP_DIR=`mktemp --tmpdir=$INIT_DIR -d`

cd $TMP_DIR

if [ "x$BUILD_R" != "x" ]; then
  # BUILD_R is tru
  curl -sSL -o tmp.tar.gz --retry 10 http://ftp.heanet.ie/mirrors/cran.r-project.org/src/base/R-3/R-3.1.3.tar.gz
  mkdir $TMP_DIR/R-build
  tar -C $TMP_DIR/R-build --strip-components 1 -zxf tmp.tar.gz
  cd $TMP_DIR/R-build
  ./configure --enable-R-shlib --with-cairo=yes --prefix=$INST_PATH
  make -j$CPU
  make check
  make install
  cd $TMP_DIR
fi

curl -sSL https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_1.0-3.tar.gz > VGAM_1.0-3.tar.gz
Rscript $INIT_DIR/libInstall.R $R_LIBS_USER

cd $INIT_DIR
rm -rf $TMP_DIR
