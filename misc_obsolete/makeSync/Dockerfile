
# packaging dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        apt-utils \
        build-essential \
        autotools-dev \
        autoconf \
        automake \
        software-properties-common \
        libncurses5-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        openjdk-8-jdk \
        maven \
        wget \
        sed \
        gzip \
        perl \
        gnupg2 \
        libopenblas-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev


# install Python
RUN apt-get install -y --no-install-recommends \
        python3-dev \
        python3-setuptools \
        python3-pip \
        python2.7-dev \
        python-setuptools \
        python-pip \
        python-pysam && \
        pip install setuptools --upgrade && \
        pip install pip --upgrade && \
        pip3 install setuptools --upgrade && \
        pip3 install pip --upgrade


# install HTSlib
ARG HTSLIB_VERSION=1.9
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz && \
    tar xzf ${HTSLIB_VERSION}.tar.gz && \
    cd htslib-${HTSLIB_VERSION} && \
    autoheader && \
    autoreconf && \
    ./configure --prefix /opt/htslib/${HTSLIB_VERSION} && \
    make && \
    make install && \
    cd .. && \
    rm ${HTSLIB_VERSION}.tar.gz && \
    rm -rf htslib-${HTSLIB_VERSION}
ENV PATH="/opt/htslib/${HTSLIB_VERSION}/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/htslib/${HTSLIB_VERSION}/lib:${LD_LIBRARY_PATH}"


# install SAMTools
ARG SAMTOOLS_VERSION=1.9
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/archive/${SAMTOOLS_VERSION}.tar.gz && \
    tar xzf ${SAMTOOLS_VERSION}.tar.gz && \
    cd samtools-${SAMTOOLS_VERSION} && \
    autoheader && \
    autoreconf -Wno-syntax && \
    ./configure \
       --prefix /opt/samtools/${SAMTOOLS_VERSION} \
       --with-htslib=/opt/htslib/${HTSLIB_VERSION} && \
    make && \
    make install && \
    cd /opt && \
    rm -rf ${SAMTOOLS_VERSION}.tar.gz && \
    rm -rf samtools-${SAMTOOLS_VERSION}
ENV PATH="/opt/samtools/${SAMTOOLS_VERSION}/bin:${PATH}"
