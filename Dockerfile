# syntax=docker/dockerfile:1
FROM continuumio/anaconda3:latest

RUN apt update
RUN apt-get install -y git wget libncurses-dev libbz2-dev liblzma-dev unzip vim bzip2
RUN apt-get install -y gcc g++ zlib1g-dev build-essential
RUN apt-get install -y openjdk-11-jdk
RUN apt-get install -y r-base

WORKDIR /

# build SAMTools
RUN wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
RUN tar xvjf samtools-1.6.tar.bz2
WORKDIR /samtools-1.6
RUN ./configure
RUN make && make install
WORKDIR /

# build GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
RUN unzip gatk-4.3.0.0.zip
WORKDIR /gatk-4.3.0.0
RUN sed -e 's/python/python3/g' < gatk > gatk4
RUN mv gatk4 gatk
RUN chmod a+x gatk

# build BCFTools, this includes htslib, tabix and bgzip
WORKDIR /
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
RUN tar xvjf bcftools-1.18.tar.bz2
WORKDIR /bcftools-1.18
RUN ./configure
RUN make && make install
WORKDIR /bcftools-1.18/htslib-1.18
RUN make && make install

# BWA 0.7.17
WORKDIR /
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /bwa
RUN make

# snpEff-5.1d, unpacks to "snpEff"
WORKDIR /
RUN wget https://networks.systemsbiology.net/downloads/snpEff-5.1d.tar.gz
RUN tar xfz snpEff-5.1d.tar.gz

# Python dependencies
RUN pip install globalsearch

# Boto3 for AWS integration
RUN conda install -y boto3

# SRA Tools for FASTQ conversion
RUN conda install -y -c conda-forge ossuuid

# just install SRA-tools directly, Anaconda gives me too many headaches
WORKDIR /
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar xvfz sratoolkit.current-ubuntu64.tar.gz
ENV PATH="$PATH:/sratoolkit.3.0.10-ubuntu64/bin"

# Sickle
RUN conda install -y -c bioconda sickle-trim

# TB Profiler:
# Installed in its own environment since it uses dependencies that
# are incompatible with the ones in here, most prominently it depends on samtools 1.12
#
# NOTE That in order to use this, we need to
# $ conda activate TBprofiler2
# $ tb-profiler
# in python.subprocess, conda activate does not without an interactive shell !!!
# So use conda run -n TBprofiler2 --live-stream tb-profiler ...
WORKDIR /
RUN wget https://networks.systemsbiology.net/downloads/tbp_v4.3.0_linux.txt
RUN conda create --name TBprofiler2 --file tbp_v4.3.0_linux.txt

# Install pipeline software and data
RUN git clone https://github.com/baliga-lab/resr_pipeline.git
WORKDIR /resr_pipeline
#RUN wget https://networks.systemsbiology.net/downloads/bwa_mtb_reference_genome-20231107.tar.gz
#RUN tar xvf bwa_mtb_reference_genome-20231107.tar.gz
#RUN rm bwa_mtb_reference_genome-20231107.tar.gz
RUN wget https://networks.systemsbiology.net/downloads/MTB_ancestor_reference-20240122.tar.gz
RUN tar xvf MTB_ancestor_reference-20240122.tar.gz
RUN rm MTB_ancestor_reference-20240122.tar.gz

# Install AWS CLI
WORKDIR /
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -u awscliv2.zip
RUN ./aws/install

# Final: start the container in the bwa_pipeline directory
WORKDIR /resr_pipeline
