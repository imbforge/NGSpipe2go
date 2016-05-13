##
## BUILD image for the RNAseq pipeline
## AVAILABLE COMMANDS (still to be run as root):
##    build: 
##      $ docker build -t imbforge/rnaseq:v1 .
##    push to the docker hub: 
##      $ docker push imbforge/rnaseq:v1
##    run a container: 
##      $ docker run --rm -v ${WORKDIR}:${WORKDIR} -w ${WORKDIR} -t imbforge/rnaseq:v1 \
##         -n ${MAX_PARALLEL_PROCS} ${WORKDIR}/rnaseq_v1.2.txt ${WORKDIR}/rawdata/*.fastq.gz
##    run an interactive shell:
##      $ docker run --entrypoint=/bin/bash -ti imbforge/rnaseq:v1 -s
##
#FROM ubuntu:14.04
FROM debian:stable
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
#RUN groupadd -r imbcf && \ 
#	useradd -r -g imbcf -s /bin/bash -c "Docker image user" imbcf

## install all deps in a single RUN statement to save space
RUN	apt-get update && apt-get dist-upgrade -y && \
	apt-get install -y apt-utils && \
	apt-get install -y build-essential gcc g++ gfortran make && \
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
COPY ./deps/R-3.2.2 /usr/local/src/R-3.2.2
WORKDIR /usr/local/src/R-3.2.2
RUN	./configure --prefix=/opt/R/3.2.2 --with-readline=no && \
	make && \
	make install && \
	/opt/R/3.2.2/bin/Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("DESeq2","edgeR","rtracklayer","GenomicAlignments","GenomicRanges","Rsubread","dupRadar"),ask=FALSE)' && \
	/opt/R/3.2.2/bin/Rscript -e 'install.packages(c("ggplot2","reshape","RColorBrewer","gplots","knitr"),dependencies=TRUE,repos="http://cran.rstudio.org")' && \
	rm -rf /usr/local/src/R-3.2.2 && \
	ln -s /opt/R/3.2.2 /opt/R/latest && \
	ln -s /opt/R/3.2.2 /opt/R/default && \
	echo "export PATH=/opt/R/3.2.2/bin:\$PATH" > /opt/R/3.2.2/env.sh && \
    chmod ugo+rx /opt/R/3.2.2/env.sh

## install FastQC
COPY ./deps/FastQC /opt/fastqc/0.11.3
RUN	ln -s /opt/fastqc/0.11.3 /opt/fastqc/latest && \
	ln -s /opt/fastqc/0.11.3 /opt/fastqc/default && \
	echo "export PATH=/opt/fastqc/0.11.3:\$PATH" > /opt/fastqc/0.11.3/env.sh && \
    chmod ugo+rx /opt/fastqc/0.11.3/env.sh

## install htseq-count
COPY ./deps/HTSeq-0.6.1p1 /usr/local/src/HTSeq-0.6.1p1
WORKDIR /usr/local/src/HTSeq-0.6.1p1
RUN mkdir -p /opt/htseq/0.6.1p1/lib/python2.7/site-packages/ && \
	export PYTHONPATH=/opt/htseq/0.6.1p1/lib/python2.7/site-packages:\$PYTHONPATH && \
	python setup.py build && \
	python setup.py install --prefix=/opt/htseq/0.6.1p1 && \
	rm -rf /usr/local/src/HTSeq-0.6.1p1 && \
	ln -s /opt/htseq/0.6.1p1 /opt/htseq/latest && \
	ln -s /opt/htseq/0.6.1p1 /opt/htseq/default && \
	echo "export PATH=/opt/htseq/0.6.1p1/bin:\$PATH" > /opt/htseq/0.6.1p1/env.sh && \
	echo "export PYTHONPATH=/opt/htseq/0.6.1p1/lib/python2.7/site-packages:\$PYTHONPATH" >> /opt/htseq/0.6.1p1/env.sh && \
    chmod ugo+rx /opt/htseq/0.6.1p1/env.sh

## install Subread
COPY ./deps/subread-1.4.6-p2 /opt/subread/1.4.6-p2
RUN ln -s /opt/subread/1.4.6-p2 /opt/subread/latest && \
	ln -s /opt/subread/1.4.6-p2 /opt/subread/default && \
	echo "export PATH=/opt/subread/1.4.6-p2/bin:\$PATH" > /opt/subread/1.4.6-p2/env.sh && \
    chmod ugo+rx /opt/subread/1.4.6-p2/env.sh

## install STAR
COPY ./deps/STAR_2.4.2a /opt/star/2.4.2a
RUN ln -s /opt/star/2.4.2a /opt/star/latest && \
	ln -s /opt/star/2.4.2a /opt/star/default && \
	echo "export PATH=/opt/star/2.4.2a/bin/Linux_x86_64_static:\$PATH" > /opt/star/2.4.2a/env.sh && \
    chmod ugo+rx /opt/star/2.4.2a/env.sh

## install samtools
COPY ./deps/samtools-1.2 /usr/local/src/samtools-1.2
WORKDIR /usr/local/src/samtools-1.2
RUN	make && \
	make install prefix=/opt/samtools/1.2 && \
	rm -rf /usr/local/src/samtools-1.2 && \
	ln -s /opt/samtools/1.2 /opt/samtools/latest && \
	ln -s /opt/samtools/1.2 /opt/samtools/default && \
	echo "export PATH=/opt/samtools/1.2/bin:\$PATH" > /opt/samtools/1.2/env.sh && \
    chmod ugo+rx /opt/samtools/1.2/env.sh

## install Picard tools
COPY ./deps/picard-tools-1.123 /opt/picard/1.123
RUN ln -s /opt/picard/1.123 /opt/picard/latest && \
	ln -s /opt/picard/1.123 /opt/picard/default && \
	echo "export PATH=/opt/picard/1.123:\$PATH" > /opt/picard/1.123/env.sh && \
    chmod ugo+rx /opt/picard/1.123/env.sh

## install RSeQC
COPY ./deps/RSeQC-2.4 /usr/local/src/RSeQC-2.4
WORKDIR /usr/local/src/RSeQC-2.4
RUN mkdir -p /opt/RSeQC/2.4/lib/python2.7/site-packages/ && \
	export PYTHONPATH=/opt/RSeQC/2.4/lib/python2.7/site-packages:\$PYTHONPATH && \
	python setup.py build && \
	python setup.py install --root=/opt/RSeQC/2.4 && \
	rm -rf /opt/RSeQC/2.4/lib && mv /opt/RSeQC/2.4/usr/local/* /opt/RSeQC/2.4 && rm -r /opt/RSeQC/2.4/usr && \
	rm -rf /usr/local/src/RSeQC-2.4 && \
	ln -s /opt/RSeQC/2.4 /opt/RSeQC/latest && \
	ln -s /opt/RSeQC/2.4 /opt/RSeQC/default && \
	echo "export PATH=/opt/RSeQC/2.4/bin:\$PATH" > /opt/RSeQC/2.4/env.sh && \
	echo "export PYTHONPATH=/opt/RSeQC/2.4/lib/python2.7/dist-packages:\$PYTHONPATH" >> /opt/RSeQC/2.4/env.sh && \
	echo "source /opt/R/3.2.2/env.sh" >> /opt/RSeQC/2.4/env.sh && \
    chmod ugo+rx /opt/RSeQC/2.4/env.sh

## install bedtools
COPY ./deps/bedtools2-2.25.0 /usr/local/src/bedtools2-2.25.0
WORKDIR /usr/local/src/bedtools2-2.25.0
RUN make && \
	make install prefix=/opt/BEDTools/2.25.0 && \
	rm -rf /usr/local/src/bedtools2-2.25.0 && \
	ln -s /opt/BEDTools/2.25.0 /opt/BEDTools/latest && \
	ln -s /opt/BEDTools/2.25.0 /opt/BEDTools/default && \
	echo "export PATH=/opt/BEDTools/2.25.0/bin:\$PATH" > /opt/BEDTools/2.25.0/env.sh && \
    chmod ugo+rx /opt/BEDTools/2.25.0/env.sh

## install UCSC tools
COPY ./deps/ucsc-2014-02-05 /opt/ucsc/2014-02-05
RUN ln -s /opt/ucsc/2014-02-05 /opt/ucsc/latest && \
	ln -s /opt/ucsc/2014-02-05 /opt/ucsc/default && \
	echo "export PATH=/opt/ucsc/2014-02-05/:\$PATH" > /opt/ucsc/2014-02-05/env.sh && \
    chmod ugo+rx /opt/ucsc/2014-02-05/env.sh

## install bpipe
COPY ./deps/bpipe-0.9.8.7 /opt/bpipe/0.9.8.7
RUN ln -s /opt/bpipe/0.9.8.7 /opt/bpipe/latest && \
	ln -s /opt/bpipe/0.9.8.7 /opt/bpipe/default && \
	echo "export PATH=\${TOOL_DEPENDENCIES}/bpipe/0.9.8.7/bin/:\$PATH" > /opt/bpipe/0.9.8.7/env.sh && \
    chmod ugo+rx /opt/bpipe/0.9.8.7/env.sh

## install pipelines and wrappers
COPY ./deps/imb-forge /opt/imb-forge
RUN ln -s /opt/imb-forge/pipelines/rnaseq/1.2 /opt/imb-forge/pipelines/rnaseq/latest && \
	ln -s /opt/imb-forge/pipelines/rnaseq/1.2 /opt/imb-forge/pipelines/rnaseq/default

## post install tasks
WORKDIR /
RUN rm -rf /usr/local/src/*
#USER imbcf
