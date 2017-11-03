FROM broadinstitute/viral-baseimage:0.1.4

LABEL maintainer "Chris Tomkins-Tinch <tomkinsc@broadinstitute.org>"

# to build:
#   docker build .
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
# Set it up so that this slow & heavy build layer is cached
# unless the requirements* files or the install scripts actually change
COPY requirements-conda.txt requirements-conda-tests.txt requirements-py3.txt $VIRAL_NGS_PATH/
COPY docker/install-viral-ngs.sh $VIRAL_NGS_PATH/docker/
COPY easy-deploy-script/easy-deploy-viral-ngs.sh $VIRAL_NGS_PATH/easy-deploy-script/
WORKDIR $INSTALL_PATH
RUN $VIRAL_NGS_PATH/docker/install-viral-ngs.sh

# Copy all of the source code into the repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_NGS_PATH/

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; source /opt/miniconda/bin/activate /opt/miniconda; echo -n 'viral-ngs version: '; $VIRAL_NGS_PATH/assembly.py --version"

# Volume setup: make external tools and data available within the container
VOLUME ["/gatk", "/novoalign", "/user-data"]
ENV GATK_PATH="/gatk" NOVOALIGN_PATH="/novoalign" VIRAL_NGS_DOCKER_DATA_PATH="/user-data"

# It's a wrapper script to load the viral-ngs environment via the easy-deploy script
# and then run any commands desired
ENTRYPOINT ["/opt/viral-ngs/source/docker/env_wrapper.sh"]
CMD ["/bin/bash"]
