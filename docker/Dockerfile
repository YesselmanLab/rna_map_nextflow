################## BASE IMAGE ######################

FROM continuumio/miniconda3:latest

################## METADATA ######################
LABEL base.image="continuumio/miniconda3:latest"
LABEL version="1"
LABEL software="rna-map-nextflow"
LABEL software.version="0.1.0"


################## MAINTAINER ######################
MAINTAINER Joe Yesselman <jyesselm@unl.edu>

################## INSTALLATION ######################

USER root

RUN /bin/bash -c "mkdir /data && mkdir /data2 && mkdir /data3 && mkdir /data4 "

## install fastqc and unzip
RUN apt-get update && \
    apt-get -y -q install build-essential zlib1g-dev fastqc gcc unzip curl openjdk-17-jdk openjdk-17-jre

RUN git clone https://github.com/najoshi/sabre.git && \
    cd sabre && \
    make && \
    cp sabre /bin

## install seqkit
RUN conda install -c bioconda seqkit

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/bin

## install bowtie2
RUN wget -q $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

## setup correct python and packages
RUN pip install cython colorlog click future plotly matplotlib pandas kaleido fastqsplitter

## install cutadapt
RUN wget -q https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip -O $DST/v1.18.zip && \
    unzip $DST/v1.18.zip -d $DST && \
    rm $DST/v1.18.zip
RUN /bin/bash -c "cd $DST/cutadapt-1.18 && python setup.py install"

## install trim_galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore $DST && \
    rm -rf TrimGalore-0.6.6 trim_galore.tar.gz


COPY python/rna_map_nextflow /bin/rna_map_nextflow 
RUN cd /bin/rna_map_nextflow && \
    pip install . 

WORKDIR /data
