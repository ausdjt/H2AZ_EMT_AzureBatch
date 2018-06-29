FROM ubuntu:16.04
RUN apt-get clean
RUN apt-get update 

# Install dependancies
RUN apt-get install -y --fix-missing --no-install-recommends \
unzip \
build-essential \
wget \
make \ 
bzip2 \
libncurses5-dev \
zlib1g-dev \ 
libbz2-dev \
software-properties-common \
debconf-utils \
python-pip \
python-numpy \
python-matplotlib \
python-pysam \
python-htseq

# Add Java support
RUN add-apt-repository ppa:webupd8team/java
RUN apt-get update 
RUN echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | debconf-set-selections
RUN apt-get install -y oracle-java8-set-default

# Folder for tools
RUN mkdir /home/apps

# Add kallisto 
RUN wget --no-check-certificate -O /home/apps/kallisto_linux-v0.44.0.tar.gz https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
RUN cd /home/apps ; tar zxvf kallisto_linux-v0.44.0.tar.gz ; mv kallisto_linux-v0.44.0 kallisto_v0.44.0

# Add picard 
RUN mkdir /home/apps/picard
RUN wget --no-check-certificate -O /home/apps/picard/picard.jar https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar

RUN apt-get install -y --fix-missing --no-install-recommends liblzma-dev

# Add Samtools
RUN wget --no-check-certificate  -O /home/apps/samtools-1.8.tar.bz2 https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
RUN cd /home/apps ; bunzip2 samtools-1.8.tar.bz2 ; tar -xvf samtools-1.8.tar ; cd samtools-1.8 ; ./configure ; make ; make install

# Add Bowtie
RUN wget --no-check-certificate  -O /home/apps/bowtie2-2.3.4.1-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip/download
RUN cd /home/apps ; unzip bowtie2-2.3.4.1-linux-x86_64.zip; 

# upgrade pip
RUN pip install --upgrade pip

# Add cutadapt
RUN apt-get install -y --fix-missing --no-install-recommends python-setuptools
RUN pip install cutadapt

# Add deeptools
RUN apt-get install -y --fix-missing --no-install-recommends git
RUN cd /home/apps ; git clone https://github.com/deeptools/deepTools ; cd /home/apps/deepTools ; python setup.py install

# Add STAR
RUN wget --no-check-certificate  -O /home/apps/2.6.0a.tar.gz https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
RUN cd /home/apps ; tar zxvf 2.6.0a.tar.gz ; cd STAR-2.6.0a/source ; make STAR


# Add Snakemake
RUN apt-get install -y --fix-missing --no-install-recommends snakemake
RUN cd /home/apps ; rm *.tar* ; rm *.zip
RUN cd /home/apps ; mv STAR-2.6.0a STAR ; mv bowtie2-2.3.4.1-linux-x86_64 bowtie2; mv kallisto_v0.44.0 kallisto ; mv samtools-1.8 samtools
RUN echo 'PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/apps"' > /etc/environment

