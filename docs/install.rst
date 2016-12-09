Installation
============


Manual Installation
-------------------


Install Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use viral-ngs, you need to install the `Conda package manager <http://conda.pydata.org/miniconda.html>`_ which is most easily obtained via the Miniconda Python distribution. Miniconda can be installed on your system without admin priviledges.

After installing Miniconda for your platform, be sure to update it::

  conda update -y conda

Configure Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viral-ngs software and its dependencies are distributed through the bioconda channel for the conda package manager. It is necessary to add this channel to the conda config::

  conda config --add channels bioconda
  conda config --add channels r
  conda config --add channels conda-forge

Make a conda environment and install viral-ngs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to install viral-ngs into its own conda environment. This ensures its dependencies do not interfere with other conda packages installed on your system. A new conda environment can be created with the following command, which will also install conda::

  conda create -n viral-ngs-env viral-ngs

Activate the viral-ngs environment and complete the install
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to finish installing viral-ngs, you will need to activate its conda environment::

  source activate viral-ngs-env

Due to license restrictions, the viral-ngs conda package cannot distribute and install GATK directly. To fully install GATK, you must download a licensed copy of GATK `from the Broad Institute <https://www.broadinstitute.org/gatk/download/>`_, and call "gatk-register," which will copy GATK into your viral-ngs conda environment::

  # (download licensed copy of GATK)
  gath-register /path/to/GenomeAnalysisTK.jar

The single-threaded version of `Novoalign <http://www.novocraft.com/products/novoalign/>`_ is installed by default. If you have a license for Novoalign to enable multi-threaded operation, viral-ngs will copy it to the viral-ngs conda environment if the ``NOVOALIGN_LICENSE_PATH`` environment variable is set. Alternatively, the conda version of Novoalign can be overridden if the ``NOVOALIGN_PATH`` environment variable is set. If you obtain a Novoalign license after viral-ngs has already been installed, it can be added to the conda environment by calling::

  # obtain a Novoalign license file: novoalign.lic
  novoalign-register-license /path/to/novoalign.lic


Activating viral-ngs once installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After viral-ngs has been installed, only one command is needed to load the environment and all of its dependencies. This is the command that must be run each time before using viral-ngs::

  source activate viral-ngs-env

To deactivate the conda environment::

  source deactivate

Easy deployment script for viral-ngs
------------------------------------

**viral-ngs** can be deployed with help of a shell script, ``easy-deploy/easy-deploy-viral-ngs.sh``. This script will install an independent copy of viral-ngs from the latest source, install all dependencies, and make it simple to activate the viral-ngs environment and create projects.  The script is available from the repository `broadinstitute/viral-ngs-deploy <https://github.com/broadinstitute/viral-ngs-deploy/tree/master/easy-deploy-script>`_.


One-line install command 
~~~~~~~~~~~~~~~~~~~~~~~~~

After downloading the easy-install shell script, this one-line command will install viral-ngs on a 64-bit macOS or Linux system::

  ./easy-deploy-script/easy-deploy-viral-ngs.sh setup

One-line install command for Broad Institute users
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This one-line command will download the ``easy-deploy-viral-ngs.sh`` script and setup viral-ngs in the current working directory. Simply ssh to one of the Broad login nodes and paste this command::

  wget https://raw.githubusercontent.com/broadinstitute/viral-ngs-deploy/master/easy-deploy-script/easy-deploy-viral-ngs.sh && chmod a+x ./easy-deploy-viral-ngs.sh && reuse UGER && qrsh -l h_vmem=10G -cwd -N "viral-ngs_deploy" -q interactive ./easy-deploy-viral-ngs.sh setup

**Note:** The script will run the install on a UGER interactive node, so you must have the ability to create to start a new interactive session. A project can be specified via ``qrsh -P "<project_name>"``

Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Installation**

* ``./easy-deploy-viral-ngs.sh setup`` Installs a fresh copy of viral-ngs, installs all dependencies, and creates a directory, ``viral-ngs-etc/``, in the current working directory.

Resulting directories::

  viral-ngs-etc/
      conda-env/
      viral-ngs/
      mc3/

**Activating the environment**

* ``source ./easy-deploy-viral-ngs.sh load`` Loads the dotkits needed by viral-ngs and activates the Python virtual environment

**Creating a project directory**

* ``./easy-deploy-viral-ngs.sh create-project <project_name>`` Creates a directory for a new Snakemake-compatible project, with data directories and symlinked run scripts. Copies in the files ``Snakefile`` and ``config.yaml``


Resulting directories::

  viral-ngs-analysis-software/
      projects/
          <project_name>/
              Snakefile
              bin/ (symlink)
              config.yaml
              data/
              log/
              reports/
              run-pipe_LSF.sh (symlink)
              run-pipe_UGER.sh (symlink)
              samples-assembly-failures.txt
              samples-assembly.txt
              samples-depletion.txt
              samples-runs.txt
              tmp/
              venv/ (symlink)
              [...other project files...]


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
``NOVOALIGN_LICENSE_PATH`` environment variables.

The easy deployment script has been tested to run on OS X 10.11 (El Capitan) and
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

Note that the EC2 instance created by the easy-deploy script is currently configured to be an m4.2xlarge, which costs ~$0.55/hour to run. It is suggested that the instance be terminated via the AWS web console once processing with viral-ngs is complete. See the `AWS page for current pricing <https://aws.amazon.com/ec2/pricing/>`_ .

Limitations
~~~~~~~~~~~

As viral-ngs does not currently build a depletion database for BMTagger or BLAST automatically,
it is the responsibility of the user to create a depletion database for use within the virtualized
viral-ngs environment. It can be created within the virtual machine (VM), or uploaded
after the fact via ``rsync``.

Running Easy Deploy
~~~~~~~~~~~~~~~~~~~

Running Easy Deploy to create a virtualized viral-ngs environment is as simple as running ``easy-deploy-virtualized/run.sh``. Before running this script, copy any data you wish to have in the vm to the ``easy-deploy-virtualized/data`` directory on your local machine. During setup, the
files will be copied into the ``~/data/`` directory of virtual machine.

To start, the script ``run.sh`` installs the necessary dependencies on the user's machine (ansible, vagrant, virtualbox, and virtualbox-aws). The provisioning is handled by Ansible, with Vagrant handling creation of the VMs and EC2 instances. On OSX it depends on Homebrew, and will install it if it is not present. It depends on having apt on linux. Ruby >=2.0 is required for vagrant-aws, so versions of Ubuntu older than 15.04 (notably 14.04 LTS) will need to have ruby >=2.0 installed and made default.

Details on Easy Deploy
~~~~~~~~~~~~~~~~~~~~~~

Per the Vagrantfile, local VM RAM usage is set to 8GB. On EC2 it currently uses an m4.2xlarge instance with 32GB of RAM and 8 vCPUs.

Ansible clones the master branch of viral-ngs from GitHub, creates a Python 3 virtual environment, and installs the viral-ngs Python dependencies. The viral-ngs tool unit tests are run to download, install, and build all of the viral-ngs tools. A ``Snakefile`` for viral-ngs is copied to the home directory of the VM (locally: ``/home/vagrant/``, on EC2: `/home/ubuntu/`), along with an associated ``config.yaml`` file. Files to contain sample names (``sample-depletion.txt``, etc.) are also created. A directory is created within the VM, ``~/data/``, to store data to be processed. This directory on the VM is synced to the ``./data/`` directory on the host machine, relative to the location of the ``easy-deploy-virtualized/Vagrantfile``. On local VMs, syncing of the directory is two-way and fast. On EC2 instances, the syncing is currently one way (local->EC2) due to Vagrant limitations.
