FROM broadinstitute/viral-baseimage:0.1

LABEL maintainer "Chris Tomkins-Tinch <tomkinsc@broadinstitute.org>"

# to build:
#   docker build --build-arg VERSION=`git describe --tags --always` --rm .
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

# This is used to lock in the version string.
#ARG VERSION

# Silence some warnings about Readline. Checkout more over here:
# https://github.com/phusion/baseimage-docker/issues/58
ENV DEBIAN_FRONTEND noninteractive

# copy basic files
COPY docker/env_wrapper.sh docker/install-skeleton.sh docker/easy-deploy-script/easy-deploy-viral-ngs.sh /opt/docker/
RUN chmod a+x /opt/docker/*.sh


##############################
# Prepare viral-ngs user and installation directory
##############################
ENV INSTALL_PATH="/opt/viral-ngs"
ENV VIRAL_NGS_PATH="/opt/viral-ngs/source"
RUN mkdir -p $VIRAL_NGS_PATH
COPY *.py tools util pipes requirements* test .git $VIRAL_NGS_PATH/
WORKDIR $INSTALL_PATH
RUN /opt/docker/install-skeleton.sh


# Volume setup: make external tools and data available within the container
VOLUME ["/gatk", "/novoalign", "/user-data"]
ENV GATK_PATH="/gatk" NOVOALIGN_PATH="/novoalign" VIRAL_NGS_DOCKER_DATA_PATH="/user-data"

# Silence some warnings about Readline. Checkout more over here:
# https://github.com/phusion/baseimage-docker/issues/58
ENV DEBIAN_FRONTEND teletype

# It's a wrapper script to load the viral-ngs environment via the easy-deploy script
# and then run any commands desired
ENTRYPOINT ["/opt/docker/env_wrapper.sh"]
CMD ["/bin/bash"]
