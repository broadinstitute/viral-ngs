Installation
============


Cloud compute implementations
-----------------------------

Docker Images
~~~~~~~~~~~~~

To facilitate cloud compute deployments, we publish a complete Docker 
image with associated dependencies to the Docker registry at `quay.io 
<https://quay.io/repository/broadinstitute/viral-ngs>`_. Simply ``docker 
pull quay.io/broadinstitute/viral-ngs`` for the latest stable version.


DNAnexus
~~~~~~~~

This assembly pipeline is also available via the DNAnexus cloud
platform. RNA paired-end reads from either HiSeq or MiSeq instruments
can be securely uploaded in FASTQ or BAM format and processed through
the pipeline using graphical and command-line interfaces. Instructions
for the cloud analysis pipeline are available at
https://github.com/dnanexus/viral-ngs/wiki


Google Cloud Platform: deploy to GCE VM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The docker image referenced above can be directly `deployed to a Google Compute Engine VM on startup <https://cloud.google.com/compute/docs/containers/deploying-containers>`_. The main things you will need to do are:

* Make sure to allocate a larger-than-default root disk for the VM. Google's Container Optimized OS defaults to a very small disk which is not large enough to unpack our Docker image. Increase to at least 20GB (or more if you want to localize data).
* When setting up the VM for launch, make sure you open the "Advanced container options" hidden options and select "Allocate a buffer for  STDIN" and "Allocate a pseudo-TTY" before launching. Otherwise you won't be able to ssh into them!
* Sometimes you will need to invoke "bash" manually upon login to get the correct environment.


Google Cloud Platform: dsub
~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the command line functions in viral-ngs are accessible from the docker image_ and can be invoked directly using dsub_.

.. _dsub: https://cloud.google.com/genomics/v1alpha2/dsub
.. _image: https://quay.io/repository/broadinstitute/viral-ngs

Here is an example invocation of ``illumina.py illumina_demux`` (replace the project with your GCP project, and the input, output-recursive, and logging parameters with URIs within your GCS buckets)::

  dsub --project my-google-project-id --zones "us-central1-*" \
    --image quay.io/broadinstitute/viral-ngs \
    --name illumina_demux \
    --logging gs://mybucket/logs \
    --input FC_TGZ=gs://mybucket/flowcells/160907_M04004_0066_000000000-AJH8U.tar.gz \
    --output-recursive OUTDIR=gs://mybucket/demux \
    --command 'illumina.py illumina_demux ${FC_TGZ} 1 ${OUTDIR}' \
    --min-ram 30 \
    --min-cores 8 \
    --disk-size 2000

The speed of disk write and read operations is linearly proportional to the disk size, hitting the maximum disk speed somewhere around 1-8TB (depending on your I/O pattern). See `GCE documentation <https://cloud.google.com/compute/docs/disks/performance>`_.


Manual Installation
-------------------


Install Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use viral-ngs, you need to install the `Conda package manager <http://conda.pydata.org/miniconda.html>`_ which is most easily obtained via the Miniconda Python distribution. Miniconda can be installed on your system without admin priviledges.

After installing Miniconda for your platform, be sure to update it::

  conda update -y conda

Configure Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viral-ngs software and its dependencies are distributed through the a channel of the conda package manager. It is necessary to add this channel to the conda config::

  conda config --add channels r
  conda config --add channels defaults 
  conda config --add channels conda-forge 
  conda config --add channels bioconda
  conda config --add channels broad-viral

Make a conda environment and install viral-ngs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to install viral-ngs into its own conda environment. This ensures its dependencies do not interfere with other conda packages installed on your system. A new conda environment can be created with the following command, which will also install conda::

  conda create -n viral-ngs-env viral-ngs

Activate the viral-ngs environment and complete the install
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to finish installing viral-ngs, you will need to activate its conda environment::

  source activate viral-ngs-env

Due to license restrictions, the viral-ngs conda package cannot distribute and install GATK directly. To fully install GATK, you must download a licensed copy of GATK v3.8 `from the Broad Institute <https://software.broadinstitute.org/gatk/download/archive>`_, and call "gatk3-register," which will copy GATK into your viral-ngs conda environment::

  mkdir -p /path/to/gatk_dir
  wget -O - 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.6-0-g89b7209' | tar -xjvC /path/to/gatk_dir
  gatk3-register /path/to/gatk_dir/GenomeAnalysisTK.jar

The single-threaded version of `Novoalign <http://www.novocraft.com/products/novoalign/>`_ is installed by default. If you have a license for Novoalign to enable multi-threaded operation, viral-ngs will copy it to the viral-ngs conda environment if the ``NOVOALIGN_LICENSE_PATH`` environment variable is set. Alternatively, the conda version of Novoalign can be overridden if the ``NOVOALIGN_PATH`` environment variable is set. If you obtain a Novoalign license after viral-ngs has already been installed, it can be added to the conda environment by calling::

  # obtain a Novoalign license file: novoalign.lic
  novoalign-license-register /path/to/novoalign.lic

Activating viral-ngs once installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After viral-ngs has been installed, only one command is needed to load the environment and all of its dependencies. This is the command that must be run each time before using viral-ngs::

  source activate viral-ngs-env

To deactivate the conda environment::

  source deactivate

Easy deployment script for viral-ngs
------------------------------------

**viral-ngs** can be deployed with help of a shell script, ``easy-deploy/easy-deploy-viral-ngs.sh``. This script will install an independent copy of viral-ngs from the latest source, install all dependencies, and make it simple to activate the viral-ngs environment and create projects.  The script is available from the repository `broadinstitute/viral-ngs <https://github.com/broadinstitute/viral-ngs/tree/master/easy-deploy-script>`_.

After downloading the script and making it executable, viral-ngs can be installed on a 64-bit macOS or Linux system via::

  ./easy-deploy-viral-ngs.sh setup

One-line install command 
~~~~~~~~~~~~~~~~~~~~~~~~~

This one-line command will install viral-ngs on a 64-bit macOS or Linux system::

  wget https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/easy-deploy-script/easy-deploy-viral-ngs.sh && chmod a+x ./easy-deploy-viral-ngs.sh && ./easy-deploy-viral-ngs.sh setup

One-line install command for Broad Institute cluster users
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This one-line command will download the ``easy-deploy-viral-ngs.sh`` script and setup viral-ngs in the current working directory. Simply ssh to one of the Broad login nodes and paste this command::

  wget https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/easy-deploy-script/easy-deploy-viral-ngs.sh && chmod a+x ./easy-deploy-viral-ngs.sh && reuse UGER && qrsh -l h_vmem=10G -cwd -N "viral-ngs_deploy" -q interactive ./easy-deploy-viral-ngs.sh setup

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

