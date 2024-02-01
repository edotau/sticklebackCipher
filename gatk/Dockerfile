FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    coreutils \
    samtools \
    bwa \
    openjdk-17-jdk \
    openjdk-17-jre \
    fastqc \
    python3 \
    python3-pip \
    && pip3 install --upgrade pip \
    && pip3 install cutadapt \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy packages for bbmap, gatk, and trim_galore
COPY ./pkg/BBMap_39.06.tar.gz /usr/local/bin/bbmap/
COPY ./pkg/gatk-4.5.0.0.tar.gz /usr/local/bin/gatk-4/
COPY ./pkg/trim_galore.tar.gz /usr/local/bin/trim_galore/

# Extract and clean up in one layer
RUN tar -xvf /usr/local/bin/bbmap/BBMap_39.06.tar.gz -C /usr/local/bin/bbmap/ \
    && rm /usr/local/bin/bbmap/BBMap_39.06.tar.gz \
    && tar -xvf /usr/local/bin/gatk-4/gatk-4.5.0.0.tar.gz -C /usr/local/bin/gatk-4/ \
    && rm /usr/local/bin/gatk-4/gatk-4.5.0.0.tar.gz \
    && tar -xvf /usr/local/bin/trim_galore/trim_galore.tar.gz -C /usr/local/bin/trim_galore/ \
    && rm /usr/local/bin/trim_galore/trim_galore.tar.gz

# Set environment variables for tool paths
ENV BBDUK_PATH=/usr/local/bin/bbmap/bbmap/ \
    TRIM_GALORE_PATH=/usr/local/bin/trim_galore/TrimGalore-0.6.10/ \
    GATK_PATH=/usr/local/bin/gatk-4/gatk-4.5.0.0/

# Copy source code or other resources
COPY ./src ./src

# Optional: Specify default command or entrypoint
CMD ["bash"]
