FROM quay.io/broadinstitute/viral-baseimage:0.2.2

LABEL maintainer "viral-ngs@broadinstitute.org"

# to build:
#   git describe --tags --always --dirty > VERSION
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
	CONDA_DEFAULT_ENV=viral-ngs-env \
	CONDA_ENVS_PATH="$MINICONDA_PATH/envs"
ENV \
	PATH="$VIRAL_NGS_PATH:$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV/bin:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
	CONDA_PREFIX=$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV \
	PYTHONPATH=$VIRAL_NGS_PATH \
	JAVA_HOME=$MINICONDA_PATH

# Prepare viral-ngs user and installation directory
# Set it up so that this slow & heavy build layer is cached
# unless the requirements* files or the install scripts actually change
WORKDIR $INSTALL_PATH
RUN conda create -n $CONDA_DEFAULT_ENV python=3.10
RUN echo "source activate $CONDA_DEFAULT_ENV" > ~/.bashrc
RUN hash -r
COPY docker $VIRAL_NGS_PATH/docker/
COPY requirements-conda.txt requirements-conda-tests.txt $VIRAL_NGS_PATH/
RUN $VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_NGS_PATH/requirements-conda.txt $VIRAL_NGS_PATH/requirements-conda-tests.txt
RUN $VIRAL_NGS_PATH/docker/install-gatk.sh

# Copy all of the source code into the repo
# This probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time.
COPY util $VIRAL_NGS_PATH/util/
COPY tools $VIRAL_NGS_PATH/tools/
COPY *.py VERSION* $VIRAL_NGS_PATH/

RUN git config --global --add safe.directory $VIRAL_NGS_PATH

# R fails unless you do this, CollectInsertSizeMetrics needs R, why is conda R broken
RUN	ln -s /lib/x86_64-linux-gnu/libreadline.so.7 /lib/x86_64-linux-gnu/libreadline.so.6

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; read_utils.py --version"

CMD ["/bin/bash"]
