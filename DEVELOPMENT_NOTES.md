## Developer Documentation
This page lists information for developers working on viral-ngs.

### Installation, dependencies, and manual deployment requirements

#### Dependency install destinations
When Python and binary dependencies for viral-ngs are installed by conda, they can end up in several locations. The default and preferred method of installation assumes a conda environment is active in the current shell, complete with [environment variables we can access to specify the path of the active environment](https://github.com/broadinstitute/viral-ngs/blob/master/tools/__init__.py#L240). In this case, conda packages are installed in the active conda environment. If conda is installed and available on the path but no environment is currently active, viral-ngs dependencies are installed in isolation within `viral-ngs/tools/build/conda-tools/{default}` (unless this location is overridden in the CondaPackage() constructor). For tools without a conda recipe (as may be the case on certain platforms, like Mac OSX), or where conda install fails, custom install methods are used to download and build some tools.

#### Adding a new tool or dependency
When adding a new tool or dependency to viral-ngs, check to see if a conda package is already available either on the default channel (`conda search <package_name>`), or on the bioconda channel (`conda search -c bioconda <package_name>`). If so, it will needed to be added to the conda recipe template for viral-ngs. If a recipe is unavailable, it will first need to be added to a particular conda channel. [Bioconda](https://github.com/bioconda/bioconda-recipes) is used by default.

#### Changing dependency versions
The viral-ngs package installed by `conda install viral-ngs` from the [broad-viral channel](https://anaconda.org/broad-viral/viral-ngs) depends on a conda build recipe distributed in this repository. The recipe files source the various Python and binary depedencies of viral-ngs as conda packages, including version numbers, from the `requirements-*.txt` files within this repository.

#### Automated deployment and release actions
When a new tagged version of viral-ngs is [released](https://github.com/broadinstitute/viral-ngs/releases), the conda package will be updated automatically by a TravisCI deploy hook in the test build for the tagged branch. First the recipe will be updated to reflect the new version number, source archive, and dependency list. Next TravisCI will build a package for the recipe via `conda build` and upload it to the [broad-viral channel](https://anaconda.org/broad-viral/viral-ngs) of anaconda.org. If the `conda build` is successful, a remote build will be triggered on TravisCI for the [broadinstitute/viral-ngs-deploy](https://github.com/broadinstitute/viral-ngs-deploy) repository in order to build and upload a Docker image to the [broadinstitute/viral-ngs](https://hub.docker.com/r/broadinstitute/viral-ngs/) repository on Docker Hub.

### (Automated) testing 
[Travis CI](https://travis-ci.org/broadinstitute/viral-ngs) performs automated unit and integration tests for viral-ngs on each branch and pull request. Unit tests are run on each new branch commit, and longer integration tests are performed on pull requests to help ensure the stability of the `master` branch. Pull requests are gated to ensure merging to `master` is allowed only if all tests pass. The Travis configuration is specified in `.travis.yml`, and relies on files stored within `viral-ngs/travis/`.

A few notes on testing: 
- Travis Linux with Py3 is the authoritative test server
- Tox is present for convenience (local use only, not currently in use on Travis)
- `py.test` is used in place of `nose` or the built-in `unittest` framework for its generative testing, fixtures, and parallelized execution of tests. 
- `.flake8`, `.pylintrc`, and `.style.yapf` are available in the repository root, and should be used when running the associated tools
- During Travis tests an encrypted tarball of third-party resources is downloaded into the build environment to provide a licensed copy of GATK and Novoalign. For security, forks of viral-ngs [will not have](https://docs.travis-ci.com/user/pull-requests#Security-Restrictions-when-testing-Pull-Requests) the encrypted resources available for testing on Travis, nor will pull requests from forks.
- Many tests for viral-ngs depend on static input files distributed with the repository. The test input files reside in `viral-ngs/test/input/<TestClassName>`. Within specific unit tests, class-specific test input files may be accessed via the function `util.file.get_test_input_path(self)`. The parent directory for all static test files can be accessed via `util.file.get_test_path()`

### Building documentation
Documentation is built automatically for certain branches of viral-ngs by [Read the Docs](http://viral-ngs.readthedocs.io/en/latest/). The documentation template files reside within `viral-ngs/docs`, and are formatted in standard docutils [reStructuredText format](http://docutils.sourceforge.net/rst.html). [Pandoc](http://pandoc.org/) may be used for converting from other formats (such as Markdown) to reStructuredText. The `sphinx-argparse` module is used to automatically generate documentation for the argparse parsers used in viral-ngs.
