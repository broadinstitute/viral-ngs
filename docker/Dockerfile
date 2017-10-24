FROM phusion/baseimage:0.9.22

LABEL maintainer "Chris Tomkins-Tinch <tomkinsc@broadinstitute.org>"

# to build:
#   docker build --rm .
# to run:
#   Download licensed copies of GATK and Novoalign to the host machine (for Linux-64)
#   export GATK_PATH=/path/to/gatk/
#   export NOVOALIGN_PATH=/path/to/novoalign/
#   docker run --rm -v $GATK_PATH:/gatk -v $NOVOALIGN_PATH:/novoalign -v /path/to/dir/on/host:/user-data -t -i <image_ID> "<command>.py subcommand"
# if you receive a "no space on device" error:
#   docker kill $(docker ps -a -q)
#   docker rm $(docker ps -a -q)
#   docker rmi $(docker images -q)
#   docker volume rm $(docker volume ls -qf dangling=true)

ARG VIRAL_NGS_VERSION

# Silence some warnings about Readline. Checkout more over here:
# https://github.com/phusion/baseimage-docker/issues/58
ENV DEBIAN_FRONTEND noninteractive

##############################
# System packages, Google Cloud SDK, and locale
##############################
# removing /var/lib/apt/lists/* frees some space
# crcmod needed for Google Cloud SDK
# 
RUN apt-get update \
    && apt-get install -y -qq --no-install-recommends ca-certificates wget rsync curl bzip2 python less nano vim locales gcc python-dev python-setuptools \
    && apt-get upgrade -y \
    && curl -sSL https://sdk.cloud.google.com | bash -s - "--install-dir=/opt &> /dev/null" \
    && easy_install -Uq pip \
    && pip install -Uq crcmod \
    && apt-get remove -y -qq gcc python-dev python-setuptools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && locale-gen en_US.UTF-8 \
    && hash -r

ENV PATH $PATH:/opt/google-cloud-sdk/bin

# Set default locale to en_US.UTF-8
ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

# grab gosu for easy step-down from root
ENV GOSU_VERSION 1.10
RUN set -ex; \
    fetchDeps=' \
        ca-certificates \
        wget \
    '; \
    apt-get update; \
    apt-get install -y --no-install-recommends $fetchDeps; \
    rm -rf /var/lib/apt/lists/*; \
    dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')"; \
    wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch"; \
    wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch.asc"; \
# verify the signature
    export GNUPGHOME="$(mktemp -d)"; \
    gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4; \
    gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu; \
    rm -r "$GNUPGHOME" /usr/local/bin/gosu.asc; \
    chmod +x /usr/local/bin/gosu; \
# verify that the binary works
    gosu nobody true;

##############################
# Prepare user and installation directory
##############################
ENV INSTALL_DIR="/opt/viral-ngs"
RUN  mkdir -p $INSTALL_DIR

##############################
# Volume setup
##############################
# make external tools and data available within the container
VOLUME ["/gatk", "/novoalign", "/user-data"]
ENV GATK_PATH="/gatk" NOVOALIGN_PATH="/novoalign" VIRAL_NGS_DOCKER_DATA_PATH="/user-data"
# It's a wrapper script to load the viral-ngs environment via the easy-deploy script
# and then run any commands desired
# Put this up here to avoid long waiting time after heavier build layer below
ADD ./env_wrapper.sh $INSTALL_DIR/env_wrapper.sh
RUN chmod -R +x $INSTALL_DIR/env_wrapper.sh

##############################
# Setup viral-ngs
##############################
# download and use the viral-ngs easy deploy script to install viral-ngs
#RUN curl -LO https://raw.githubusercontent.com/broadinstitute/viral-ngs-deploy/master/easy-deploy-script/easy-deploy-viral-ngs.sh
# Docker does not support conditional copies, but we can get around it by copying a file we know to exist and then wildcarded files
COPY Dockerfile viral-ngs*.tar.bz2 $INSTALL_DIR/
ADD ./easy-deploy-script/easy-deploy-viral-ngs.sh $INSTALL_DIR/easy-deploy-viral-ngs.sh
RUN chmod +x $INSTALL_DIR/easy-deploy-viral-ngs.sh
WORKDIR $INSTALL_DIR

ENV VIRAL_NGS_VERSION $VIRAL_NGS_VERSION
# if a *.tar.bz2 file exists in the INSTALL_DIR, assume it is a viral-ngs package and install it
RUN echo "Running easy install with viral-ngs version $VIRAL_NGS_VERSION $(ls -1 ./*.tar.bz2 2> /dev/null | head -1)"
RUN bash -c 'echo $(if [ "$VIRAL_NGS_VERSION" == "" ]; then echo "./easy-deploy-viral-ngs.sh setup `ls -1 ./*.tar.bz2 2> /dev/null | head -1`"; else echo "./easy-deploy-viral-ngs.sh setup --viral-ngs-version $VIRAL_NGS_VERSION"; fi)'
RUN bash -c '$(if [ "$VIRAL_NGS_VERSION" == "" ]; then echo "./easy-deploy-viral-ngs.sh setup `ls -1 ./*.tar.bz2 2> /dev/null | head -1`"; else echo "./easy-deploy-viral-ngs.sh setup --viral-ngs-version $VIRAL_NGS_VERSION"; fi)'

ENV DEBIAN_FRONTEND teletype

ENTRYPOINT ["/opt/viral-ngs/env_wrapper.sh"]
CMD ["/bin/bash"]
