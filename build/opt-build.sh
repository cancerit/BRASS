#! /bin/bash

set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
set -u

## vcftools
if [ ! -e $SETUP_DIR/vcftools.success ]; then
  curl -sSL --retry 10 -o distro.tar.gz https://github.com/vcftools/vcftools/releases/download/v${VER_VCFTOOLS}/vcftools-${VER_VCFTOOLS}.tar.gz
  rm -rf distro/*
  tar --strip-components 2 -C distro -xzf distro.tar.gz
  cd distro
  ./configure --prefix=$INST_PATH --with-pmdir=lib/perl5
  make -j$CPU
  make install
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/vcftools.success
fi

### cgpVcf
if [ ! -e $SETUP_DIR/cgpVcf.success ]; then
  curl -sSL --retry 10 -o distro.tar.gz https://github.com/cancerit/cgpVcf/archive/${VER_CGPVCF}.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/cgpVcf.success
fi

## bedtools
if [ ! -e $SETUP_DIR/bedtools.success ]; then
  curl -sSL --retry 10 -o $INST_PATH/bin/bedtools https://github.com/arq5x/bedtools2/releases/download/v${VER_BEDTOOLS}/bedtools.static.binary
  chmod +x $INST_PATH/bin/bedtools
  chmod -w $INST_PATH/bin/bedtools
  touch $SETUP_DIR/bedtools.success
fi


## blat
if [ ! -e $SETUP_DIR/blat.success ]; then
  curl -sSL --retry 10 -o blat.zip https://hgwdev.gi.ucsc.edu/~kent/src/blatSrc${VER_BLAT}.zip
  unzip -qu blat.zip
  cd $SETUP_DIR/blatSrc
  BINDIR=$SETUP_DIR/blat/bin
  export BINDIR
  export MACHTYPE
  mkdir -p $BINDIR
  make -j$CPU
  cp $BINDIR/blat $INST_PATH/bin/.
  cd $SETUP_DIR
  rm -rf blat.zip blatSrc
  touch $SETUP_DIR/blat.success
fi

## ssearch36
if [ ! -e $SETUP_DIR/fasta36.success ]; then
  curl -sSL --retry 10 -o distro.tar.gz https://github.com/wrpearson/fasta36/releases/download/fasta-v${VER_FASTA36}/fasta-${VER_FASTA36}-linux64.tar.gz
  tar -zxf distro.tar.gz ./fasta-${VER_FASTA36}/bin/ssearch36
  cp ./fasta-${VER_FASTA36}/bin/ssearch36 $OPT/bin/.
  rm -rf distro.tar.gz
  touch $SETUP_DIR/fasta36.success
fi

### VAGrENT
if [ ! -e $SETUP_DIR/VAGrENT.success ]; then
  curl -sSL --retry 10 -o distro.tar.gz https://github.com/cancerit/VAGrENT/archive/${VER_VAGRENT}.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/VAGrENT.success
fi

### grass
if [ ! -e $SETUP_DIR/grass.success ]; then
  curl -sSL --retry 10 -o distro.tar.gz https://github.com/cancerit/grass/archive/${VER_GRASS}.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro

  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/grass.success
fi
