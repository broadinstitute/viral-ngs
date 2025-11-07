## Developer Documentation
This page lists information for developers working on viral-ngs.

### Modifying code and testing

The current dev, build, and deploy paradigm is intentionally docker-centric. This
means that developers will need docker working within their dev environment--and
not much else (other than git, and a text/code editor). The code base is also
modularized and layered. In order to work on code changes, you will:

1. check out the git code repository for this module on your local host machine (`git clone https://github.com/broadinstitute/viral-core.git`) and edit with your favorite code/text editor
1. docker `pull` and `run` the image `FROM` which this is built, while mounting your local git checkout into the container (`docker run -it --rm -v $(pwd):/opt/viral-ngs/source quay.io/broadinstitute/viral-core`)
1. if your dev branch has altered the conda dependencies in `requirements-conda.txt`, make sure to `conda install` them (more detailed instrucitons pending). You may (optionally) want to snapshot the new docker image locally if you want to continue using it and skip this step in the future (`docker commit <image hash> local/viral-core-dev`)
1. test code and execution interactively within the container (`cd /opt/viral-ngs/source; pytest -rsxX -n auto test/unit`)
1. push changes back to github (from your host machine) for automated CI testing & builds, using standard, collaborative github code review processes

### Machinery under the hood

#### Dependency install destinations
When Python and binary dependencies for viral-ngs are installed by conda, they can end up in several locations. The default and preferred method of installation assumes a conda environment is active in the current shell, complete with [environment variables we can access to specify the path of the active environment](https://github.com/broadinstitute/viral-ngs/blob/master/tools/__init__.py#L240). In this case, conda packages are installed in the active conda environment. 

#### Adding a new tool or dependency
When adding a new tool or dependency to viral-ngs, check to see if a conda package is already available either on the default channel (`conda search <package_name>`), or on the bioconda channel (`conda search -c bioconda <package_name>`). If so, it will needed to be added to the conda recipe template for viral-ngs. If a recipe is unavailable, it will first need to be added to a particular conda channel. [Bioconda](https://github.com/bioconda/bioconda-recipes) is used by default.


### (Automated) testing 
GitHub Actions performs automated unit and integration tests for viral-ngs on each branch and pull request. Unit tests are run on each new branch commit, and longer integration tests are performed on pull requests to help ensure the stability of the `master` branch. Pull requests are gated to ensure merging to `master` is allowed only if all tests pass. The test configuration is specified in `.github/workflows/build.yml`, and relies on files stored within `viral-ngs/github_actions_ci/`.

A few notes on testing: 
- The docker container is the authoritative test, dev, and production execution environment.
- `py.test` is used in place of `nose` or the built-in `unittest` framework for its generative testing, fixtures, and parallelized execution of tests. 
- Many tests for viral-ngs depend on static input files distributed with the repository. The test input files reside in `viral-ngs/test/input/<TestClassName>`. Within specific unit tests, class-specific test input files may be accessed via the function `util.file.get_test_input_path(self)`. The parent directory for all static test files can be accessed via `util.file.get_test_path()`

#### The build matrix
Each commit on any branch, and any pull request, will trigger a build. Branch commits will test code from a specific commit hash. Pull requests will test the simulated result of merging a branch HEAD onto the target branch HEAD. For each build, the following jobs are launched:
1. Docker container first. A docker image is built and deployed to the Docker registry at quay.io. Master branch images are pushed to `quay.io/broadinstitute/viral-core:latest` and are also given a versioned tag. Non-master branch images and pull requests are pushed to `quay.io/broadinstitute/viral-core` with versioned tags. The docker build is preceded by a docker pull of the docker image associated with the previous build parental to this commit in order to utilize layer caching. Note that our tool dependencies result in a very large docker image (1GB compressed, this is about 10x the typical size for a docker image). The Dockerfile builds the tool dependencies before incorporating the full viral-core source code. This means that most docker image builds will be extremely fast: usually 10-20 seconds. The docker push/deploy is similarly fast, since the Docker registry already has most of the layers, and only the new source code layer needs to upload. The docker pull of the 1GB image takes about 2 minutes, so altogether this step takes about ~6 minutes. However, if your code commit alters anything in `requirements-*.txt` or the easy deploy script, it will rebuild the heavy conda install layer, adding another 10 minutes or so to this build. The docker push requires login credentials for a docker registry (e.g. DockerHub, Quay.io, GCP, AWS), stored as variables in the GitHub actions repository secrets.
1. Documentation is built automatically. It is not deployed to Read the Docs from GitHub Actions; the Actions test only exists in order to bring the developer's attention to any auto-build problems. Read the Docs has its own auto-build process separate from GitHub Actions (see section below) but it does not notify anyone of its build failures. This usually completes in less than 1 minute.
1. The docker container is pulled and `py.test` is run inside the container. Integration and unit tests are run with every branch commit. Coverage reports are copied from container to host and sent to coveralls.io from the Actions host VM.

### Building documentation
Documentation is built automatically for master branches of viral-core by [Read the Docs](http://viral-core.readthedocs.io/en/latest/). The documentation template files reside within `viral-core/docs`, and are formatted in standard docutils [reStructuredText format](http://docutils.sourceforge.net/rst.html). [Pandoc](http://pandoc.org/) may be used for converting from other formats (such as Markdown) to reStructuredText. The `sphinx-argparse` module is used to automatically generate documentation for the argparse parsers used in viral-ngs.
