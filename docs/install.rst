Installation
============


Manual Installation
-------------------

System dependencies
~~~~~~~~~~~~~~~~~~~

This is known to install cleanly on most modern Linux systems with Python,
Java, and some basic development libraries.  On Ubuntu 14.04 LTS, the
following APT packages should be installed on top of the vanilla setup::

  python3 python3-pip python3-nose
  python-software-properties

.. (comment out below)
  zlib zlib1g zlib1g-dev
  libblas3gf libblas-dev liblapack3gf liblapack-dev
  libatlas-dev libatlas3-base libatlas3gf-base libatlas-base-dev
  gfortran
  oracle-java8-installer
  libncurses5-dev

.. (comment out below)
 The Fortran libraries (including blas and atlas) are required to install
 numpy via pip from source. numpy is not actually required if you have
 Python 3.4, if you want to avoid this system dependency.

**Java >= 1.7** is required by GATK and Picard.


Python dependencies
~~~~~~~~~~~~~~~~~~~

The **command line tools** require Python >= 2.7 or >= 3.4. Required packages
(pysam and Biopython) are listed in requirements.txt and can be
installed the usual pip way::

  pip install -r requirements.txt

Additionally, in order to use the **pipeline infrastructure**, Python 3.4
is required (Python 2 is not supported) and you must install snakemake
as well::

  pip install -r requirements-pipes.txt

However, most of the real functionality is encapsulated in the command line
tools, which can be used without any of the pipeline infrastructure.

You should either sudo pip install or use a virtualenv (recommended).


Tool dependencies
~~~~~~~~~~~~~~~~~

A lot of effort has gone into writing auto download/compile wrappers for
most of the bioinformatic tools we rely on here. 

Most tools will attemp a conda-based install first, before falling back to an install handled entirely by our wrappers. To make use of the conda-based install, you will need to have Anaconda or miniconda installed on your system:

http://conda.pydata.org/docs/install/quick.html#miniconda-quick-install-requirements

The tools will auto-download and install the first time they are needed by any command. If you want to pre-install all of the external tools, simply type this::

  nosetests -v test.unit.test_tools

However, there are two tools in particular that cannot be auto-installed
due to licensing restrictions.  You will need to download and install
these tools on your own (paying for it if your use case requires it) and
set environment variables pointing to their installed location.

 * GATK - http://www.broadinstitute.org/gatk/
 * Novoalign - http://www.novocraft.com/products/novoalign/

The environment variables you will need to set are ``GATK_PATH`` and
``NOVOALIGN_PATH``. These should be set to the full directory path
that contains these tools (the jar file for GATK and the executable
binaries for Novoalign).

In order to run GATK, you will need to have an appropriate version of 
the Java JDK installed. As of this writing, Java 1.7 is required for 
GATK 3.3.0. 

Alternatively, if you are using the Snakemake pipelines, you can create
a dictionary called "env_vars" in the config.yaml file for Snakemake,
and the pipelines will automatically set all environment variables prior
to running any scripts.

The version of MOSAIK we use seems to fail compile on GCC-4.9 but compiles
fine on GCC-4.4. We have not tried intermediate versions of GCC, nor the
latest versions of MOSAIK.

Virtualized Installation (Easy Deploy)
--------------------------------------

The viral-ngs package includes a script that can be used to set up a complete virtualized 
environment for running viral-ngs either on a local machine via VirtualBox, or on AWS EC2. 
This is an easiesr way to get the software up and running, as it sets up most 
dependencies automatically within an environment known to work.

Requirements
~~~~~~~~~~~~

As noted above, GATK and NovoAlign cannot be installed automatically due to 
licensing restrictions. In order to run the easy deployment script, you will
first need to license and download these tools, and set the ``GATK_PATH`` and 
``NOVOALIGN_PATH`` environment variables. 

The easy deployment script has been tested to run on OS X 10.10 (Yosemite) and
Ubuntu 15.04 (Vivid Vervet).


Requirements for running on AWS EC2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to deploy a virtualized viral-ngs environment to AWS EC2, you will first need
to set up the appropriate credentials for creating EC2 instances. AWS credentials and 
SSH keypairs are passed in as environment variables, and ``run.sh`` will prompt for 
the values if the environment variables are not set (though the values given 
interactively are ephemeral).

The following environment variables are needed:

 * ``EC2_ACCESS_KEY_ID``
 * ``EC2_SECRET_ACCESS_KEY``
 * ``EC2_REGION`` (ex. "us-west-2")
 * ``EC2_KEYPAIR_NAME`` (ex. "my-ssh-keypair")
 * ``EC2_PRIVATE_KEY_PATH`` (ex. "my-ssh-keypair.pem")
 * ``EC2_SECURITY_GROUP`` (ex. "ssh-only-group")

For more information, see the following AWS pages:

* `Getting set up with AWS <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html>`_
* `How to create an AWS EC2 key pair <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair>`_
* `Defining security group rules <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-network-security.html#adding-security-group-rule>`_
* `List of EC2 regions <https://docs.aws.amazon.com/general/latest/gr/rande.html#ec2_region>`_
 
Note that the EC2 instance created by the easy-deploy script is currently configured to be an m4.2xlarge, which costs ~$0.55/hour to run. It is suggested that the instance be terminated via the AWS web console once processing with viral-ngs is comple. See the `AWS page for current pricing <https://aws.amazon.com/ec2/pricing/>`_ .

Limitations
~~~~~~~~~~~

As viral-ngs does not currently build a depletion database for BMTagger or BLAST automatically, 
it is the responsibility of the user to create a depletion database for use within the virtualized
viral-ngs environment. It can be created within the virtual machine (VM), or uploaded
after the fact via ``rsync``.

Running Easy Deploy
~~~~~~~~~~~~~~~~~~~

Running Easy Deploy to create a virtualized viral-ngs environment is as simple as running ``easy-deploy/run.sh``. Before running this script, copy any data you wish to have in the vm to the ``easy-deploy/data`` directory on your local machine. During setup, the 
files will be copied into the ``~/data/`` directory of virtual machine.

To start, the script ``run.sh`` installs the necessary dependencies on the user's machine (ansible, vagrant, virtualbox, and virtualbox-aws). The provisioning is handled by Ansible, with Vagrant handling creation of the VMs and EC2 instances. On OSX it depends on Homebrew, and will install it if it is not present. It depends on having apt on linux. Ruby >=2.0 is required for vagrant-aws, so versions of Ubuntu older than 15.04 (notably 14.04 LTS) will need to have ruby >=2.0 installed and made default.

Details on Easy Deploy
~~~~~~~~~~~~~~~~~~~~~~

Per the Vagrantfile, local VM RAM usage is set to 8GB. On EC2 it currently uses an m4.2xlarge instance with 32GB of RAM and 8 vCPUs.

Ansible clones the master branch of viral-ngs from GitHub, creates a Python 3 virtual environment, and installs the viral-ngs Python dependencies. The viral-ngs tool unit tests are run to download, install, and build all of the viral-ngs tools. A ``Snakefile`` for viral-ngs is copied to the home directory of the VM (locally: ``/home/vagrant/``, on EC2: `/home/ubuntu/`), along with an associated ``config.yaml`` file. Files to contain sample names (``sample-depletion.txt``, etc.) are also created. A directory is created within the VM, ``~/data/``, to store data to be processed. This directory on the VM is synced to the ``./data/`` directory on the host machine, relative to the location of the ``easy-deploy/Vagrantfile``. On local VMs, syncing of the directory is two-way and fast. On EC2 instances, the syncing is currently one way (local->EC2) due to Vagrant limitations.