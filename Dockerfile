FROM broadinstitute/viral-baseimage:0.1.2

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

# DEBIAN_FRONTEND: Silence some warnings about Readline. Read more over here:
# https://github.com/phusion/baseimage-docker/issues/58
ENV DEBIAN_FRONTEND=noninteractive INSTALL_PATH="/opt/viral-ngs" VIRAL_NGS_PATH="/opt/viral-ngs/source"

# copy basic files
COPY docker/env_wrapper.sh docker/install-viral-ngs.sh easy-deploy-script/easy-deploy-viral-ngs.sh $INSTALL_PATH/
RUN chmod a+x $INSTALL_PATH/*.sh

# Prepare viral-ngs user and installation directory
COPY . $VIRAL_NGS_PATH/
WORKDIR $INSTALL_PATH
RUN ./install-viral-ngs.sh

# Volume setup: make external tools and data available within the container
VOLUME ["/gatk", "/novoalign", "/user-data"]
# DEBIAN_FRONTEND: Silence some warnings about Readline. Read more over here:
# https://github.com/phusion/baseimage-docker/issues/58
ENV GATK_PATH="/gatk" NOVOALIGN_PATH="/novoalign" VIRAL_NGS_DOCKER_DATA_PATH="/user-data" DEBIAN_FRONTEND=teletype

# It's a wrapper script to load the viral-ngs environment via the easy-deploy script
# and then run any commands desired
ENTRYPOINT ["/opt/viral-ngs/env_wrapper.sh"]
CMD ["/bin/bash"]
