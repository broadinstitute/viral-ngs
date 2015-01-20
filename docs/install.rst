Installation
============


System dependencies
-------------------

This is known to install cleanly on most modern Linux systems with Python,
Java, and some basic development libraries.  On Ubuntu 14.04 LTS, the
following APT packages should be installed on top of the vanilla setup::

  python3 python3-pip python3-nose
  python-software-properties
  zlib zlib1g zlib1g-dev
  libblas3gf libblas-dev liblapack3gf liblapack-dev
  libatlas-dev libatlas3-base libatlas3gf-base libatlas-base-dev
  gfortran
  oracle-java8-installer
  libncurses5-dev

The Fortran libraries (including blas and atlas) are required to install
numpy via pip from source. numpy is not actually required if you have
Python 3.4, if you want to avoid this system dependency.

**Java >= 1.7** is required by GATK and Picard.


Python dependencies
-------------------

The **command line tools** require Python >= 2.7 or >= 3.4. Required packages
(like pysam and Biopython) are listed in requirements.txt and can be
installed the usual pip way::

  pip install -r requirements.txt

Additionally, in order to use the **pipeline infrastructure**, Python 3.4
is required (Python 2 is not supported) and you must install snakemake
as well::

  pip install snakemake==3.2 yappi=0.94

You should either sudo pip install or use a virtualenv (recommended).


Tool dependencies
-----------------

A lot of effort has gone into writing auto download/compile wrappers for
most of the bioinformatic tools we rely on here. They will auto-download
and install the first time they are needed by any command. If you want
to pre-install all of the external tools, simply type this::

  python -m unittest test.test_tools.TestToolsInstallation -v

However, there are two tools in particular that cannot be auto-installed
due to licensing restrictions.  You will need to download and install
these tools on your own (paying for it if your use case requires it) and
set environment variables pointing to their installed location.

 * GATK - http://www.broadinstitute.org/gatk/
 * Novoalign - http://www.novocraft.com/products/novoalign/

The environment variables you will need to set are GATK_PATH and
NOVOALIGN_PATH. These should be set to the full directory path
that contains these tools (the jar file for GATK and the executable
binaries for Novoalign).

Alternatively, if you are using the Snakemake pipelines, you can create
a dictionary called "env_vars" in the config.json file for Snakemake,
and the pipelines will automatically set all environment variables prior
to running any scripts.

The version of MOSAIK we use seems to fail compile on GCC-4.9 but compiles
fine on GCC-4.4. We have not tried intermediate versions of GCC, nor the
latest versions of MOSAIK.
