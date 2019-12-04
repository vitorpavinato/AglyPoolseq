FROM broadinstitute/gatk:4.1.2.0

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


# install R
ARG DEBIAN_FRONTEND=noninteractive
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
#    apt install -y --no-install-recommends \
#        apt-transport-https software-properties-common \
RUN apt install -y --no-install-recommends \
        r-base \
        r-recommended && \

    rm -rf /var/lib/apt/lists/*


# install bam2fastq
ARG BAM2FASTQ_VERSION=1.1.0
WORKDIR /opt
RUN wget https://gsl.hudsonalpha.org/static/software/bam2fastq-${BAM2FASTQ_VERSION}.tgz && \
    tar xzf bam2fastq-${BAM2FASTQ_VERSION}.tgz && \
    mkdir bam2fastq && \
    mv bam2fastq-${BAM2FASTQ_VERSION} bam2fastq/${BAM2FASTQ_VERSION} && \
    cd bam2fastq/${BAM2FASTQ_VERSION} && \
    make && \
    cd ../.. && \
    rm bam2fastq-${BAM2FASTQ_VERSION}.tgz
ENV PATH="/opt/bam2fastq/${BAM2FASTQ_VERSION}:${PATH}"


# install bwa
ARG BWA_VERSION=0.7.15
WORKDIR /opt
RUN echo $BWA_VERSION && \
    mkdir bwa && \
    cd bwa && \
    wget https://github.com/lh3/bwa/archive/v${BWA_VERSION}.tar.gz && \
    tar xzf v${BWA_VERSION}.tar.gz && \
    mv bwa-${BWA_VERSION} ${BWA_VERSION} && \
    cd ${BWA_VERSION} && \
    make && \
    cd .. && \
    rm v${BWA_VERSION}.tar.gz
ENV PATH="/opt/bwa/${BWA_VERSION}:${PATH}"


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


# install Cutadapt
ARG CUTADAPT_VERSION=2.3
RUN python3 -m pip install cutadapt==${CUTADAPT_VERSION}

ENV JRE_HOME=${JAVA_HOME}/jre

# install Picard Tool
ARG PICARD_VERSION=2.20.2
WORKDIR /opt
RUN mkdir -p /opt/picard/${PICARD_VERSION} && \
    cd /opt/picard/${PICARD_VERSION} && \
    wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar
ENV CLASSPATH="/opt/picard/${PICARD_VERSION}/picard.jar:$CLASSPATH"
ENV PICARD="/opt/picard/${PICARD_VERSION}/picard.jar"


# install GenomeAnalysisTK
ARG GATK_VERSION=3.8-1
WORKDIR /opt/gatk/
