FROM  quay.io/wtsicgp/pcap-core:5.0.5 as builder

USER  root

# ALL tool versions used by opt-build.sh
ENV VER_BEDTOOLS="2.28.0"
ENV VER_VCFTOOLS="0.1.16"
ENV VER_BLAT="35"
ENV VER_CGPVCF="v2.2.1"
# for ssearch36, don't include 'v'
ENV VER_FASTA36="36.3.8g"
ENV VER_VAGRENT="v3.5.0"
ENV VER_GRASS="v2.1.1"

RUN apt-get -yq update

RUN apt-get install -qy --no-install-recommends lsb-release

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu `lsb_release -cs`/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get -yq update

# no benefit of combined in builder stage
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends g++
RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends gcc
RUN apt-get install -yq --no-install-recommends pkg-config
RUN apt-get install -yq --no-install-recommends zlib1g-dev
RUN apt-get install -yq --no-install-recommends unzip
RUN apt-get install -yq --no-install-recommends libpng12-dev
RUN apt-get install -yq --no-install-recommends r-base
RUN apt-get install -yq --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -yq --no-install-recommends libxml2-dev
RUN apt-get install -yq --no-install-recommends libgit2-dev
RUN apt-get install -yq --no-install-recommends libssl-dev
RUN apt-get install -yq --no-install-recommends libcairo2-dev
RUN apt-get install -yq --no-install-recommends gfortran
RUN apt-get install -yq --no-install-recommends libblas-dev
RUN apt-get install -yq --no-install-recommends libboost-all-dev
RUN apt-get install -yq --no-install-recommends libpstreams-dev

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS

# build tools from other repos
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY Rsupport Rsupport
COPY distros distros
ADD build/opt-build-local-deps.sh build/
RUN bash build/opt-build-local-deps.sh $OPT

COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:16.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="v6.3.3" \
      description="BRASS docker"

RUN apt-get -yq update \
&& apt-get install -qy --no-install-recommends lsb-release

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu `lsb_release -cs`/" >> /etc/apt/sources.list \
&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& apt-get -yq update

RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
ca-certificates \
libperlio-gzip-perl \
bzip2 \
psmisc \
time \
zlib1g \
liblzma5 \
libncurses5 \
p11-kit \
exonerate \
libcairo2 \
gfortran \
r-base \
libboost-iostreams-dev \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
