## Easy deployment of viral-ngs

**viral-ngs** can be deployed with help from the script in this directory, `easy-deploy-viral-ngs.sh`. This script will install an independent copy of viral-ngs from the latest source, install all dependencies, and make it simple to activate the viral-ngs environment and create projects.

### One-line install command 

This one-line command will install viral-ngs on a 64-bit macOS or Linux system:

```
./easy-deploy-script/easy-deploy-viral-ngs.sh setup
```

### One-line install command for Broad Institute users

This one-line command will download the `easy-deploy-viral-ngs.sh` script and setup viral-ngs in the current working directory. Simply ssh to one of the login nodes and paste this command:

    wget https://raw.githubusercontent.com/broadinstitute/viral-ngs-deploy/master/easy-deploy-script/easy-deploy-viral-ngs.sh && chmod a+x ./easy-deploy-viral-ngs.sh && reuse UGER && qrsh -l h_vmem=4G -cwd -N "viral-ngs_deploy" -q interactive ./easy-deploy-viral-ngs.sh setup

**Note:** The script will run the install on a UGER interactive node, so you must have the ability to create to start a new interactive session. A project can be specified via `qrsh -P "<project_name>"` 

### Usage

* `./easy-deploy-viral-ngs.sh setup` Installs a fresh copy of viral-ngs,  installs all dependencies, and creates a directory, `viral-ngs-etc`, in the current working directory. 

**Resulting directories**:

```
viral-ngs-etc/
    conda-env/
    viral-ngs/
    mc3/
```

* `source ./easy-deploy-viral-ngs.sh load` Loads the dotkits needed by viral-ngs and activates the Python virtual environment

* `./easy-deploy-viral-ngs.sh create-project <project_name>` Creates a directory for a new Snakemake-compatible project, with data directories and symlinked run scripts. Copies in the files `Snakefile` and `config.yaml`

**Resulting directories**:

```
viral-ngs-etc/
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
