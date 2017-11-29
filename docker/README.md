## viral-ngs Docker container

### Introduction
The parent directory contains a `Dockerfile` suitable for building a containerized instance of viral-ngs capable of running on any platform that supports [Docker](https://www.docker.com/). 

Detailed documentation for Docker is [available online](https://docs.docker.com/).

### To build
  From the viral-ngs directory, run:

  `docker build .`

  See `.travis.yml` for an example invocation of this.
  
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


### Running as a non-root user within the container

An optional entrypoint script is provided that invokes gosu to step down from root:

```docker run --rm -ti --entrypoint=/opt/viral-ngs/source/docker/gosu-entrypoint.sh <image_ID>```


### Debugging

#### Shell usage
To use a shell within a viral-ngs Docker container, run in interactive mode:

```docker run --rm -t -i <image_ID>```

#### Clean slate
If you receive a "no space on device" error, sometimes a fresh start can be helpful. You can run these commands to remove **ALL** current docker images, containers, and volumes (be careful! the commands will also remove Docker items unrelated to viral-ngs):
```shell
docker kill $(docker ps -a -q)
docker rm $(docker ps -a -q)
docker rmi $(docker images -a -q)
docker volume rm $(docker volume ls -qf dangling=true)
```
