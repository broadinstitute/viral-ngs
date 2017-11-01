FROM broadinstitute/viral-baseimage:0.1.4

LABEL maintainer "Chris Tomkins-Tinch <tomkinsc@broadinstitute.org>"

# to build:
#   docker build --rm .
#
# to run:
#   docker run --rm <image_ID> "<command>.py subcommand"
#
# to run interactively:
#   docker run --rm -it <image_ID> bash
#
# to run with GATK and/or Novoalign:
#   Download licensed copies of GATK and Novoalign to the host machine (for Linux-64)
#   export GATK_PATH=/path/to/gatk/
#   export NOVOALIGN_PATH=/path/to/novoalign/
#   docker run --rm -v $GATK_PATH:/gatk -v $NOVOALIGN_PATH:/novoalign -v /path/to/dir/on/host:/user-data <image_ID> "<command>.py subcommand"

ENV INSTALL_PATH="/opt/viral-ngs" VIRAL_NGS_PATH="/opt/viral-ngs/source"

# Prepare viral-ngs user and installation directory
COPY . $VIRAL_NGS_PATH/
WORKDIR $INSTALL_PATH
RUN $VIRAL_NGS_PATH/docker/install-viral-ngs.sh

# Volume setup: make external tools and data available within the container
VOLUME ["/gatk", "/novoalign", "/user-data"]
ENV GATK_PATH="/gatk" NOVOALIGN_PATH="/novoalign" VIRAL_NGS_DOCKER_DATA_PATH="/user-data"

# It's a wrapper script to load the viral-ngs environment via the easy-deploy script
# and then run any commands desired
ENTRYPOINT ["/opt/viral-ngs/source/docker/env_wrapper.sh"]
CMD ["/bin/bash"]
