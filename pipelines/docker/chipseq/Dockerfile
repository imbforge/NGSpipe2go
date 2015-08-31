##
## BUILD image for the ChIPseq pipeline
## AVAILABLE COMMANDS:
##    build: 
##      $ docker build -t imbforge/chipseq:v1 .
##    push to the docker hub: 
##      $ docker push imbforge/chipseq:v1
##    run a container: 
##      $ docker run --rm -v ${WORKDIR}:${WORKDIR} -w ${WORKDIR} -t imbforge/chipseq:v1
##        ${WORKDIR}/bpipe.txt ${WORKDIR}/rawdata/*.fastq.gz
##    run an interactive shell:
##      $ docker run -ti imbforge/chipseq:v1 bash
##
#FROM ubuntu:14.04
FROM debian:wheezy
MAINTAINER Sergi Sayols <s.sayolspuig@imb-mainz.de>
ENTRYPOINT ["/opt/bpipe/default/bin/bpipe","run"]

## env vars
ENV HOME /home
ENV DEBIAN_FRONTEND noninteractive
ENV TOOL_DEPENDENCIES /opt

## create user to avoid running the thing as root
# Note: Users and groups in an image get a non-deterministic UID/GID in that the
#   "next" UID/GID gets assigned regardless of image rebuilds. So, if it’s critical,
#   you should assign an explicit UID/GID.
RUN groupadd -r imbcf && \ 
	useradd -r -g imbcf -s /bin/bash -c "Docker image user" imbcf

## install all deps in a single RUN statement to save space
RUN	apt-get update && apt-get dist-upgrade -y && \
	apt-get install -y apt-utils && \
	apt-get install -y build-essential gcc g++ gfortran make && \
	apt-get install -y libboost-dev && \
	apt-get install -y xorg-dev libglu1-mesa-dev libpango1.0-dev && \ 
	apt-get install -y libxml2-dev libcurl4-openssl-dev && \
	apt-get install -y openjdk-7-jre-headless && \
	apt-get install -y python2.7-dev python-numpy python-matplotlib python-setuptools && \
	apt-get install -y zlib1g-dev libncurses5-dev && \
	apt-get install -y bc && \
	apt-get clean autoclean && \
	apt-get autoremove -y && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## install R and the needed packages
COPY ./deps/R-3.1.2 /usr/local/src/R-3.1.2
WORKDIR /usr/local/src/R-3.1.2
RUN	./configure --prefix=/opt/R/3.1.2 --with-readline=no && \
	make && \
	make install && \
	/opt/R/3.1.2/bin/Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("GenomicAlignments","GenomicRanges","ShortRead"),ask=FALSE)' && \
	/opt/R/3.1.2/bin/Rscript -e 'install.packages(c("snow","knitr"),dependencies=TRUE,repos="http://cran.rstudio.org")' && \
	rm -rf /usr/local/src/R-3.1.2 && \
	ln -s /opt/R/3.1.2 /opt/R/latest && \
	ln -s /opt/R/3.1.2 /opt/R/default && \
	echo "export PATH=/opt/R/3.1.2/bin:\$PATH" > /opt/R/3.1.2/env.sh

## install SPP
COPY ./deps/spp_1.11.tar.gz /usr/local/src/
WORKDIR /usr/local/src
RUN	/opt/R/3.1.2/bin/R CMD INSTALL spp_1.11.tar.gz && \
	rm spp_1.11.tar.gz

## install FastQC
COPY ./deps/FastQC /opt/fastqc/0.11.2
RUN	ln -s /opt/fastqc/0.11.2 /opt/fastqc/latest && \
	ln -s /opt/fastqc/0.11.2 /opt/fastqc/default && \
	echo "export PATH=/opt/fastqc/0.11.2:\$PATH" > /opt/fastqc/0.11.2/env.sh

## install bowtie
COPY ./deps/bowtie-1.1.1 /opt/bowtie/1.1.1
RUN	ln -s /opt/bowtie/1.1.1 /opt/bowtie/latest && \
	ln -s /opt/bowtie/1.1.1 /opt/bowtie/default && \
	echo "export PATH=/opt/bowtie/1.1.1:\$PATH" > /opt/bowtie/1.1.1/env.sh

## install samtools
COPY ./deps/samtools-1.2 /usr/local/src/samtools-1.2
WORKDIR /usr/local/src/samtools-1.2
RUN	make && \
	make install prefix=/opt/samtools/1.2 && \
	rm -rf /usr/local/src/samtools-1.2 && \
	ln -s /opt/samtools/1.2 /opt/samtools/latest && \
	ln -s /opt/samtools/1.2 /opt/samtools/default && \
	echo "export PATH=/opt/samtools/1.2/bin:\$PATH" > /opt/samtools/1.2/env.sh

## install Picard tools
COPY ./deps/picard-tools-1.123 /opt/picard/1.123
RUN ln -s /opt/picard/1.123 /opt/picard/latest && \
	ln -s /opt/picard/1.123 /opt/picard/default && \
	echo "export PATH=/opt/picard/1.123:\$PATH" > /opt/picard/1.123/env.sh

## install bedtools
COPY ./deps/bedtools2 /usr/local/src/bedtools2
WORKDIR /usr/local/src/bedtools2
RUN make && \
	make install prefix=/opt/BEDTools/2.23.0 && \
	rm -rf /usr/local/src/bedtools2 && \
	ln -s /opt/BEDTools/2.23.0 /opt/BEDTools/latest && \
	ln -s /opt/BEDTools/2.23.0 /opt/BEDTools/default && \
	echo "export PATH=/opt/BEDTools/2.23.0/bin:\$PATH" > /opt/BEDTools/2.23.0/env.sh

## install UCSC tools
COPY ./deps/ucsc-2014-02-05 /opt/ucsc/2014-02-05
RUN ln -s /opt/ucsc/2014-02-05 /opt/ucsc/latest && \
	ln -s /opt/ucsc/2014-02-05 /opt/ucsc/default && \
	echo "export PATH=/opt/ucsc/2014-02-05/:\$PATH" > /opt/ucsc/2014-02-05/env.sh

## install macs2
COPY ./deps/macs2-2.1.0.20150420 /usr/local/src/macs2-2.1.0.20150420
WORKDIR /usr/local/src/macs2-2.1.0.20150420
RUN	mkdir -p /opt/macs2/2.1.0.20150420/lib/python2.7/site-packages && \
	export PYTHONPATH=/opt/macs2/2.1.0.20150420/lib/python2.7/site-packages:\$PYTHONPATH && \
	python setup.py build && \
	python setup.py install --prefix=/opt/macs2/2.1.0.20150420 && \
	rm -rf /usr/local/src/macs2-2.1.0.20150420 && \
	ln -s /opt/macs2/2.1.0.20150420 /opt/macs2/latest && \
	ln -s /opt/macs2/2.1.0.20150420 /opt/macs2/default && \
	echo "export PATH=/opt/macs2/2.1.0.20150420/bin:\$PATH" > /opt/macs2/2.1.0.20150420/env.sh && \
	echo "export PYTHONPATH=/opt/macs2/2.1.0.20150420/lib/python2.7/site-packages:\$PYTHONPATH" >> /opt/macs2/2.1.0.20150420/env.sh

## install bpipe
COPY ./deps/bpipe-0.9.8.5 /opt/bpipe/0.9.8.5
RUN ln -s /opt/bpipe/0.9.8.5 /opt/bpipe/latest && \
	ln -s /opt/bpipe/0.9.8.5 /opt/bpipe/default && \
	echo "export PATH=\${TOOL_DEPENDENCIES}/bpipe/0.9.8.5/bin/:\$PATH" > /opt/bpipe/0.9.8.5/env.sh

## install pipelines and wrappers
COPY ./deps/imb-forge /opt/imb-forge
RUN ln -s /opt/imb-forge/pipelines/chipseq/1.0 /opt/imb-forge/pipelines/chipseq/latest && \
	ln -s /opt/imb-forge/pipelines/chipseq/1.0 /opt/imb-forge/pipelines/chipseq/default

## post install tasks
WORKDIR /
RUN rm -rf /usr/local/src/*
USER imbcf