FROM ubuntu:latest

RUN DEBIAN_FRONTEND=noninteractive \
    apt-get update \
    && apt-get install -y python3 python3-pip \
    && pip3 install --upgrade pip \
    && pip3 install cutadapt \
    && apt-get install -y \
    samtools bwa \
    openjdk-17-jdk \
    openjdk-17-jre \
    fastqc \
    vim \
    && rm -rf /var/lib/apt/lists/*

COPY ./pkg ./pkg

RUN mkdir /usr/local/bin/bbmap && mkdir /usr/local/bin/gatk-4 && mkdir /usr/local/bin/trim_galore

RUN tar -xvf ./pkg/BBMap_39.06.tar.gz -C /usr/local/bin/bbmap/ \
    && rm -rf ./pkg/BBMap_39.06.tar.gz

RUN tar -xvf ./pkg/gatk-4.5.0.0.tar.gz -C /usr/local/bin/gatk-4/ \
    && rm -rf ./pkg/gatk-4.5.0.0.tar.gz

RUN tar -xvf ./pkg/trim_galore.tar.gz -C /usr/local/bin/trim_galore/ \
    && rm -rf ./pkg/trim_galore.tar.gz

ENV BBDUK=/usr/local/bin/bbmap/bbmap/bbduk.sh
ENV TRIM_GALORE=/usr/local/bin/trim_galore/TrimGalore-0.6.10/trim_galore
ENV GATK=/usr/local/bin/gatk-4/gatk-4.5.0.0/gatk

COPY ./src ./src