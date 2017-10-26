## viral-ngs Docker container

### Introduction
This directory contains a `Dockerfile` suitable for building a containerized instance of viral-ngs capable of running on any platform that supports [Docker](https://www.docker.com/). 

Detailed documentation for Docker is [available online](https://docs.docker.com/).

### Considerations
The `Dockerfile` relies on the [easy-install script](https://github.com/broadinstitute/viral-ngs/tree/master/easy-deploy-script), which means that any time a new Docker image is built, it will include the latest version of viral-ngs available via the broad-viral channel of the [conda package manager](http://conda.pydata.org/docs/install/quick.html). 

In order to use all features of the viral-ngs pipeline, it is necessary to manually [license and download GATK](https://software.broadinstitute.org/gatk/). To make a licensed copy of GATK available to the Docker container, set the `$GATK_PATH` environment variable on the host machine to the path containing an extracted copy of `GenomeAnalysisTK.jar`.

By default, the pipeline will install a single-threaded version of [Novoalign](http://www.novocraft.com/products/novoalign/) from bioconda. It is possible to use a multi-threaded version to decrease compute time. To make a multi-threaded version of novoalign available to the Docker container, set the `$NOVOALIGN_PATH` environment variable on the host machine to the path containing the extracted Novocraft binaries, such as `novoalign`, `novoindex`, as well as the `novoalign.lic` license file. Note that the version of Novoalign used must support the 64-bit Linux platform of the viral-ngs Docker image, regardless of the platform hosting the Docker daemon.

### To build
  Navigate to the directory containing the `Dockerfile`, then run:

  `tar -czh . | docker build --rm -`

  The `tar` is necessary because Docker cannot dereference symlinks, and by tarring the directory, symlinks
  to files in higher filesystem levels can be used. In particular, it is assumed a symlink exists within the directory containing the `Dockerfile` to the `easy-deploy-script/` directory of this repo. Furthermore, it is a good idea to use `-t` argument to give the build a tag so it is easier to manage images later. For example:

  `tar -czh . | docker build --rm -t local/viral-ngs -`

  To build with a specific version of viral-ngs, a Docker `--build-arg` may optionally be specified. Ex.:
  ```
  tar -czh . | docker build --rm -t local/viral-ngs:1.16.0 --build-arg VIRAL_NGS_VERSION=1.16.0  -
  ```
  Note that the version of viral-ngs specified must exist in one of the channels specified in the easy-install script. As of March 2017, check available version at [broad-viral](https://anaconda.org/broad-viral/viral-ngs/files) channel.  Otherwise, build may fail with an error like this:

  ```shell
 PackageNotFoundError: Package not found: '' Package missing in current linux-64 channels:
  - viral-ngs 1.16.0*

You can search for packages on anaconda.org with

    anaconda search -t conda viral-ngs

You may need to install the anaconda-client command line client with

    conda install anaconda-client
  ```
  
### To run
Download licensed copies of GATK and Novoalign to the host machine (for Linux-64), and run:
```shell
export GATK_PATH=/path/to/gatk/
export NOVOALIGN_PATH=/path/to/novoalign/
docker run --rm -v $NOVOALIGN_PATH:/novoalign -v $GATK_PATH:/gatk -v /path/to/dir/on/host:/user-data -i <image_ID> "<command>.py subcommand"
```
The helper script `rundocker.sh` is also available to conveniently run more customized docker run commands. It requires, however, to edit the file to match with locations of GATK, Novoalign, and `viral-ngs` version. This is an example to run with `rundocker.sh`

`./rundocker.sh /path/to/datdir/on/host "<command>.py subcommand"`

Alternatively, GATK_PATH can be passed in as an argument to commands requiring it, though the jar file must be accessible within the container via a shared volume from the host. Ex:
```shell
# assumes GenomeAnalysisTK.jar is in your cwd
docker run --rm -v $(pwd):/data -i <image_ID> "assembly.py refine_assembly --GATK_PATH /data/GenomeAnalysisTK.jar"
```

The Novoalign argument, `-v $NOVOALIGN_PATH:/novoalign`, can be omitted if single-threaded Novoalign is adequate.

An arbitrary directory on the host system can be made available within the Docker container by setting a volume to mount a host directory to `/user-data` within the container. For example, to expose `~/Desktop`, add the argument: `-v ~/Desktop:/user-data`. Within the Docker container, this directory will be accessible via `/user-data`. For convenience, `/user-data` is also available via `~/data/` within the container. Other volumes from the host can be shared with the container with additional `-v` arguments.

The `<command>.py subcommand` specified must match one of the [documented viral-ngs commands](https://viral-ngs.readthedocs.io/en/latest/cmdline.html).

### Debugging

#### Shell usage
To use a shell within a viral-ngs Docker container, pass `/bin/bash` to the run command:

```docker run --rm -t -i <image_ID> /bin/bash```

Note that with the wrapper script normally used as the entrypoint, for true root access the following can be used for shell debugging of the container:

```docker run --rm -ti --entrypoint=/bin/bash <image_ID> -s```

#### Clean slate
If you receive a "no space on device" error, sometimes a fresh start can be helpful. You can run these commands to remove **ALL** current docker images, containers, and volumes (be careful! the commands will also remove Docker items unrelated to viral-ngs):
```shell
docker kill $(docker ps -a -q)
docker rm $(docker ps -a -q)
docker rmi $(docker images -a -q)
docker volume rm $(docker volume ls -qf dangling=true)
```
