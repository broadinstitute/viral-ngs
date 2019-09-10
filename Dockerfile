FROM quay.io/broadinstitute/viral-baseimage:0.1.15

LABEL maintainer "viral-ngs@broadinstitute.org"

# to build:
#   docker build .
#
# to run:
#   docker run --rm <image_ID> "<command>.py subcommand"
#
# to run interactively:
#   docker run --rm -it <image_ID>

ENV \
	INSTALL_PATH="/opt/viral-ngs" \
	VIRAL_NGS_PATH="/opt/viral-ngs/source" \
	MINICONDA_PATH="/opt/miniconda" \
	CONDA_DEFAULT_ENV=viral-ngs-env
ENV \
	PATH="$VIRAL_NGS_PATH:$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV/bin:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
	CONDA_PREFIX=$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV \
	JAVA_HOME=$MINICONDA_PATH

# Prepare viral-ngs user and installation directory
# Set it up so that this slow & heavy build layer is cached
# unless the requirements* files or the install scripts actually change
WORKDIR $INSTALL_PATH
RUN conda create -n $CONDA_DEFAULT_ENV python=3.6
RUN echo "source activate $CONDA_DEFAULT_ENV" > ~/.bashrc
RUN hash -r
COPY docker/install-viral-ngs.sh $VIRAL_NGS_PATH/docker/
COPY requirements-conda.txt requirements-conda-tests.txt requirements-py3.txt $VIRAL_NGS_PATH/
RUN $VIRAL_NGS_PATH/docker/install-viral-ngs.sh

# Copy all of the source code into the repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_NGS_PATH/

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; read_utils.py --version"

CMD ["/bin/bash"]
