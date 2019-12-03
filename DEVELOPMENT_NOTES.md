## Developer Documentation
This page lists information for developers working on viral-ngs.

### Installation, dependencies, and manual deployment requirements

#### Dependency install destinations
When Python and binary dependencies for viral-ngs are installed by conda, they can end up in several locations. The default and preferred method of installation assumes a conda environment is active in the current shell, complete with [environment variables we can access to specify the path of the active environment](https://github.com/broadinstitute/viral-ngs/blob/master/tools/__init__.py#L240). In this case, conda packages are installed in the active conda environment. If conda is installed and available on the path but no environment is currently active, viral-ngs dependencies are installed in isolation within `viral-ngs/tools/build/conda-tools/{default}` (unless this location is overridden in the CondaPackage() constructor). For tools without a conda recipe (as may be the case on certain platforms, like Mac OSX), or where conda install fails, custom install methods are used to download and build some tools.

#### Adding a new tool or dependency
When adding a new tool or dependency to viral-ngs, check to see if a conda package is already available either on the default channel (`conda search <package_name>`), or on the bioconda channel (`conda search -c bioconda <package_name>`). If so, it will needed to be added to the conda recipe template for viral-ngs. If a recipe is unavailable, it will first need to be added to a particular conda channel. [Bioconda](https://github.com/bioconda/bioconda-recipes) is used by default.

#### Changing dependency versions
The viral-ngs package installed by `conda install viral-ngs` from the [broad-viral channel](https://anaconda.org/broad-viral/viral-ngs) depends on a conda build recipe distributed in this repository. The recipe files source the various Python and binary depedencies of viral-ngs as conda packages, including version numbers, from the `requirements-*.txt` files within this repository.  When updating a package version in requirements-conda.txt, update it also in requirements-minimal.txt if it appears there as well.

#### Upgrading GATK
When upgrading the GATK to a new version:
- in requirements-conda.txt change the gatk version
- in pipes/config.yaml change the GATK_PATH to point to the correct GATK directory containing GenomeAnalysisTK.jar.
  May need to untar a .tar.bz2 from /humgen/gsa-hpprojects/GATK/bin/ into /idi/sabeti-scratch/shared-resources/software/gatk/ .
- in tools/gatk.py change TOOL_VERSION_TUPLE at the top
- in travis/install-gatk.sh change GATK_VERSION at the top
- in easy-deploy-script/easy-deploy-viral-ngs.sh 
- in docker/rundocker.sh 

### (Automated) testing 
[Travis CI](https://travis-ci.com/broadinstitute/viral-ngs) performs automated unit and integration tests for viral-ngs on each branch and pull request. Unit tests are run on each new branch commit, and longer integration tests are performed on pull requests to help ensure the stability of the `master` branch. Pull requests are gated to ensure merging to `master` is allowed only if all tests pass. The Travis configuration is specified in `.travis.yml`, and relies on files stored within `viral-ngs/travis/`.

A few notes on testing: 
- Travis Linux with Py3 is the authoritative test server
- Tox is present for convenience (local use only, not currently in use on Travis)
- `py.test` is used in place of `nose` or the built-in `unittest` framework for its generative testing, fixtures, and parallelized execution of tests. 
- `.flake8`, `.pylintrc`, and `.style.yapf` are available in the repository root, and should be used when running the associated tools
- During Travis tests an encrypted tarball of third-party resources is downloaded into the build environment to provide a licensed copy of GATK and Novoalign. For security, forks of viral-ngs [will not have](https://docs.travis-ci.com/user/pull-requests#Security-Restrictions-when-testing-Pull-Requests) the encrypted resources available for testing on Travis, nor will pull requests from forks.
- Many tests for viral-ngs depend on static input files distributed with the repository. The test input files reside in `viral-ngs/test/input/<TestClassName>`. Within specific unit tests, class-specific test input files may be accessed via the function `util.file.get_test_input_path(self)`. The parent directory for all static test files can be accessed via `util.file.get_test_path()`

#### The Travis build matrix
Each commit on any branch, and any pull request, will trigger a build on Travis CI. Branch commits will test code from a specific commit hash. Pull requests will test the simulated result of merging a branch HEAD onto the target branch HEAD. For each build, the following Travis jobs are launched:
1. Docker & WDL
   1. A docker image is built and deployed to the Docker registry at quay.io. Master branch images are pushed to `quay.io/broadinstitute/viral-ngs:latest` and are also given a versioned tag. Non-master branch images and pull requests are pushed to `quay.io/broadinstitute/viral-ngs-dev` with versioned tags. The docker build is preceded by a docker pull of the docker image associated with the previous Travis build parental to this commit in order to utilize layer caching. Note that our tool dependencies result in a very large docker image (2GB compressed, this is about 10x the typical size for a docker image). The Dockerfile builds the tool dependencies before incorporating the full viral-ngs source code. This means that most docker image builds will be extremely fast: usually 10-20 seconds. The docker push/deploy is similarly fast, since the Docker registry already has most of the layers, and only the new source code layer needs to upload. The docker pull of the 2GB image takes about 5 minutes, so altogether this step takes about 6 minutes on Travis. However, if your code commit alters anything in `requirements-*.txt` or the easy deploy script, it will rebuild the heavy conda install layer, adding another 10 minutes or so to this build. The docker push requires login credentials for a docker registry (e.g. DockerHub, Quay.io, GCP, AWS), stored as an encrypted Travis variable.
   2. After the docker image is deployed, WDL pipeline files are edited to reflect the version tag of the recently pushed docker image. A WDL validator is then run (using wdltool.jar) to ensure that all WDL files are still valid. This completes in seconds.
   3. WDL pipelines are compiled to DNAnexus workflows using dxWDL.jar. These are deployed to a DNAnexus CI project using an API token stored as an encrypted Travis variable. This completes in under a minute.
   4. A couple DNAnexus workflows are test executed in the CI project.
   4. WDL pipelines are executed with test data using Cromwell on the local Travis instance. This is a bit slow (roughly 5 mins for a simple test).
1. Documentation is built automatically. It is not deployed to Read the Docs--this test only exists on Travis in order to bring the developer's attention to any auto build problems. Read the Docs has its own auto build process separate from Travis (see section below) but it does not notify anyone of its build failures. This usually completes in less than 1 minute.
1. The `viral-ngs` conda package is built and deployed to the `broad-viral` channel. This requires anaconda.org credentials stored as an encrypted Travis variable. This takes about 10 minutes.
1. `py.test` is run on Python 3.6. Tool dependencies are installed prior to unit tests via conda. Integration and unit tests are run with every branch commit--note that this is the reverse order of the Py27 tests (unit then integration) so that errors are likely to be detected earlier in the overall build process, if they exist. The Travis cache is cleared for each tagged release, invoking a full re-install of dependencies. Normally, this job completes in 15+ minutes, about half of which is the loading of conda tool dependencies from the cache. Coverage reports are sent to coveralls.io from this Travis job only.

Some TO DO improvements for the future:
 - DNAnexus workflow testing should check output for correctness.
 - Cromwell workflow testing should check output for correctness.
 - Utilize Travis build stages.
   - All of the sub-steps of the first Docker & WDL Travis job should be broken out as separate jobs that wait for the Docker build and deploy.
   - Unit tests for Python 3.6, and possibly the conda package build, should occur within the Docker container.
   - Second-stage jobs that pull the docker image should utilize quay.io's torrent squashed image pull to reduce the time spent pulling our Docker image (currently about 5 minutes to pull from DockerHub).
   - Alternatively, we can explore creating a minimal docker image that installs only the conda pip packages (and perhaps extremely common conda tools like samtools and Picard) and leaves the rest of the conda tools out, letting them dynamically install themselves as needed using our dynamic tool install code.

### Building documentation
Documentation is built automatically for certain branches of viral-ngs by [Read the Docs](http://viral-ngs.readthedocs.io/en/latest/). The documentation template files reside within `viral-ngs/docs`, and are formatted in standard docutils [reStructuredText format](http://docutils.sourceforge.net/rst.html). [Pandoc](http://pandoc.org/) may be used for converting from other formats (such as Markdown) to reStructuredText. The `sphinx-argparse` module is used to automatically generate documentation for the argparse parsers used in viral-ngs.
