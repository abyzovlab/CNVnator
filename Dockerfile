FROM ubuntu:18.04

RUN apt update && apt install -y --no-install-recommends \
  ca-certificates \
  g++ \
  libbz2-dev \
  libcurl3-dev \
  libfreetype6-dev \
  liblzma-dev \
  libncurses5-dev \
  libreadline-dev \
  libpython2.7 \
  libz-dev \
  make \
  python-matplotlib \
  python-scipy \
  python-tk \
  curl \
 && rm -rf /var/lib/apt/lists/*

RUN curl https://root.cern/download/root_v6.18.04.Linux-ubuntu18-x86_64-gcc7.4.tar.gz | \
 tar -C /opt -xzf -

ENV PYTHONPATH=/opt/root/lib

RUN echo '/opt/root/lib' > /etc/ld.so.conf.d/root.conf \
 && ldconfig

RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | \
 tar -C /tmp -xjf - \
 && cd /tmp/samtools-* \
 && make \
 && (cd htslib-* && make mostlyclean) \
 && make mostlyclean \
 && find /tmp -name test -type d -exec rm -rf {} +

COPY ./ /tmp/CNVnator

RUN cd /tmp/CNVnator \
 && ln -s /tmp/samtools-* samtools \
 && ROOTSYS=/opt/root make \
 && mv cnvnator *.py *.pl /usr/local/bin \
 && mv pytools /usr/local/lib/python*/dist-packages \
 && cd - \
 && rm -rf /tmp/*
