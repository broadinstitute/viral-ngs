FROM ubuntu:15.04
MAINTAINER Chris Tomkins-Tinch <tomkinsc@broadinstitute.org>

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

# System packages 
RUN apt-get update && apt-get upgrade
RUN apt-get install -y curl bzip2 build-essential python wget 
RUN apt-get install -y less nano

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# add a new user
ENV DOCKER_USER viral-ngs-user
RUN useradd -ms /bin/bash $DOCKER_USER
USER $DOCKER_USER
WORKDIR /home/$DOCKER_USER

# shared volumes to make external GATK available within the container
VOLUME ["/gatk"]
ENV GATK_PATH /gatk

# shared volumes to make external Novoalign available within the container
VOLUME ["/novoalign"]
ENV NOVOALIGN_PATH /novoalign

# shared volumes to make external Novoalign available within the container
VOLUME ["/user-data"]
ENV VIRAL_NGS_DOCKER_DATA_PATH /user-data

RUN ln -s /user-data /home/$DOCKER_USER/data

# download and use the viral-ngs easy deploy script to install viral-ngs
RUN curl -LO https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/easy-deploy-script/easy-deploy-viral-ngs.sh
RUN chmod a+x easy-deploy-viral-ngs.sh
RUN bash easy-deploy-viral-ngs.sh setup

# source the viral-ngs conda environment in the .bashrc file so all
# commands are always available for a bash seshion    
RUN echo "source easy-deploy-viral-ngs.sh load &> /dev/null" >> /home/$DOCKER_USER/.bashrc

# make the entrypoint bash with login (source .bashrc)
ENTRYPOINT ["/bin/bash", "--login", "-i", "-c"]
