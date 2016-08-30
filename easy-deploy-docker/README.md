## viral-ngs Docker container

### Introduction
This directory contains a `Dockerfile` suitable for building a containerized instance of viral-ngs capable of running on any platform that supports [Docker](https://www.docker.com/). 

Detailed documentation for Docker is [available online](https://docs.docker.com/).

### Considerations
The `Dockerfile` relies on the [easy-install script](https://github.com/broadinstitute/viral-ngs/tree/master/easy-deploy-script), which means that any time a new Docker image is built, it will include the latest version of viral-ngs available via the [bioconda](https://bioconda.github.io/recipes/viral-ngs/README.html) channel of the [conda package manager](http://conda.pydata.org/docs/install/quick.html). 

In order to use all features of the viral-ngs pipeline, it is necessary to manually [license and download GATK](https://software.broadinstitute.org/gatk/). To make a licensed copy of GATK available to the Docker container, set the `$GATK_PATH` environment variable on the host machine to the path containing an extracted copy of `GenomeAnalysisTK.jar`.

By default, the pipeline will install a single-threaded version of [Novoalign](http://www.novocraft.com/products/novoalign/) from bioconda. It is possible to use a multi-threaded version to decrease compute time. To make a multi-threaded version of novoalign available to the Docker container, set the `$NOVOALIGN_PATH` environment variable on the host machine to the path containing the extracted Novocraft binaries, such as `novoalign`, `novoindex`, as well as the `novoalign.lic` license file. Note that the version of Novoalign used must support the 64-bit Linux platform of the viral-ngs Docker image, regardless of the platform hosting the Docker daemon.

### To build
  Navigate to the directory containing the `Dockerfile`, then run:
  `docker build --rm .`
  
### To run
Download licensed copies of GATK and Novoalign to the host machine (for Linux-64), and run:
```shell
export GATK_PATH=/path/to/gatk/
export NOVOALIGN_PATH=/path/to/novoalign/
docker run --rm -v $NOVOALIGN_PATH:/novoalign -v $GATK_PATH:/gatk -v /path/to/dir/on/host:/user-data -t -i <image_ID> "<command>.py subcommand"
```

The Novoalign argument, `-v $NOVOALIGN_PATH:/novoalign`, can be omitted if single-threaded Novoalign is adequate.

An arbitrary directory on the host system can be made available within the Docker container by setting a volume to mount a host directory to `/user-data` within the container. For example, to expose `~/Desktop`, add the argument: `-v ~/Desktop:/user-data`. Within the Docker container, this directory will be accessible via `/user-data`. For convenience, `/user-data` is also available via `~/data/` within the container. Other volumes from the host can be shared with the container with additional `-v` arguments.

The `<command>.py subcommand` specified must match one of the [documented viral-ngs commands](https://viral-ngs.readthedocs.io/en/latest/cmdline.html).

### Debugging

#### Shell usage
To use a shell within a viral-ngs Docker container, pass `/bin/bash` to the run command:
`docker run --rm -v $GATK_PATH:/gatk -v $NOVOALIGN_PATH:/novoalign -t -i <image_ID> /bin/bash`

#### Clean slate
If you receive a "no space on device" error, sometimes a fresh start can be helpful. You can run these commands to remove **ALL** current docker images, containers, and volumes (be careful! the commands will also remove Docker items unrelated to viral-ngs):
```shell
docker kill $(docker ps -a -q)
docker rm $(docker ps -a -q)
docker rmi $(docker images -a -q)
docker volume rm $(docker volume ls -qf dangling=true)
```
