## Easy deployment of viral-ngs for Broad Institute users

**viral-ngs** can be deployed on the Broad Institute cluster with help from the script in this directory, `viral-ngs-broad.sh`. This script will install an independent copy of viral-ngs from the latest source, install all dependencies, and make it simple to activate the viral-ngs environment and create projects.

### Dependencies

The script, `viral-ngs-broad.sh`, is intended to run on the Broad Institute cluster. It depends on the dotkits present on the Broad cluster and will not function properly in a different environment.

### One-line command

This one-line command will download the `viral-ngs-broad.sh` script and setup viral-ngs in the current working directory. Simply ssh to one of the login nodes and paste this command:

    wget https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/easy-deploy-broad/viral-ngs-broad.sh && chmod a+x ./viral-ngs-broad.sh && reuse UGER && qrsh -cwd -N "viral-ngs_deploy" -q interactive ./viral-ngs-broad.sh setup

**Note:** The script will run the install on a UGER interactive node, so you must have the ability to create to start a new interactive session. A project can be specified via `qrsh -P "<project_name>"` 

### Usage

* `viral-ngs-broad.sh setup` Installs a fresh copy of viral-ngs,  installs all dependencies, and creates a directory, `viral-ngs-analysis-software`, in the current working directory. 

**Resulting directories**:

```
viral-ngs-analysis-software/
    venv/
    viral-ngs/
```

* `source viral-ngs-broad.sh load` Loads the dotkits needed by viral-ngs and activates the Python virtual environment

* `viral-ngs-broad.sh create-project <project_name>` Creates a directory for a new Snakemake-compatible project, with data directories and symlinked run scripts. Copies in the files `Snakefile` and `config.yaml`

**Resulting directories**:

```
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
```