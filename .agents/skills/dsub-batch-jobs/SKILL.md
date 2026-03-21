# Running Batch Jobs on GCP via dsub

Use dsub to run one-off compute jobs on Google Cloud when your analysis requires
more compute/memory than the local environment, or needs specific Docker images
that are impractical to run locally.

## When to Use

- Analysis tools need >8GB RAM (e.g., VADR, BLAST, assembly)
- Need to run many independent jobs in parallel (batch processing)
- Need a specific Docker image with pre-installed tools
- Data lives in GCS and is most efficiently processed in-cloud

## Prerequisites

- **dsub** installed (ask the user where their dsub installation or venv is located)
- **gcloud CLI** authenticated with a GCP project that has Batch API enabled
- **GCS bucket** accessible by the project's default service account

## Generic Invocation

```bash
dsub --provider google-cls-v2 \
  --project <gcp-project> \
  --regions <region> \
  --image <docker-image> \
  --machine-type <machine-type> \
  --script <script.sh> \
  --tasks <tasks.tsv> \
  --logging gs://<bucket>/logs/<job-name>/
```

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--provider google-cls-v2` | Use Google Cloud Batch (recommended over Life Sciences) |
| `--project` | GCP project with Batch API enabled |
| `--regions` | Compute region (e.g., `us-central1`) |
| `--image` | Docker image to run (e.g., `staphb/vadr:1.6.4`) |
| `--machine-type` | VM type (e.g., `n1-highmem-4` for 26GB RAM) |
| `--script` | Local shell script to execute inside the container |
| `--tasks` | TSV file defining one row per job (batch mode) |
| `--logging` | GCS path for stdout/stderr logs |

## Task TSV Format

The tasks TSV defines inputs, outputs, and environment variables for each job.
Header row uses column prefixes to declare types:

```
--env VAR1	--env VAR2	--input FASTA	--output RESULT	--output LOG
value1	value2	gs://bucket/input.fasta	gs://bucket/output.txt	gs://bucket/log.txt
```

- `--env NAME` -- environment variable passed to the script
- `--input NAME` -- GCS file downloaded to a local path; the env var is set to the local path
- `--output NAME` -- local path; after the script finishes, the file is uploaded to GCS

Each non-header row is one job. All jobs run independently in parallel.

## GCP Project and Bucket Scoping

The service account running dsub jobs must have read/write access to all GCS paths
referenced in the tasks TSV. The simplest approach:

1. Use a GCP project whose default service account already has access to your data
2. Use a bucket within that same project for staging intermediate/output files
3. For ephemeral results, use a temp bucket with a lifecycle policy (e.g., 30-day auto-delete)

### Broad Viral Genomics Defaults

Most developers on the viral-ngs codebase use:
- **GCP project**: `gcid-viral-seq`
- **Staging bucket**: `gs://viral-temp-30d` (30-day auto-delete lifecycle)

These are not universal -- always confirm with the user before using them.

## Monitoring Jobs

```bash
# Check job status
dsub --provider google-cls-v2 --project <project> --jobs <job-id> --status

# Or use dstat
dstat --provider google-cls-v2 --project <project> --jobs <job-id> --status '*'

# View logs
gcloud storage cat gs://<bucket>/logs/<job-name>/<task-id>.log
```

## Tips

- **Batch over single jobs**: Always prefer `--tasks` with a TSV over individual
  `dsub` invocations. One TSV row per job is cleaner and easier to track.
- **Machine sizing**: Check your tool's memory requirements. VADR needs ~16GB;
  use `n1-highmem-4` (26GB). Most tools work fine with `n1-standard-4` (15GB).
- **Script portability**: Write the `--script` to be self-contained. It receives
  inputs/outputs as environment variables. Don't assume any local state.
- **Logging**: Always set `--logging` to a GCS path so you can debug failures.
- **Idempotency**: If re-running, dsub will create new jobs. Check for existing
  outputs before re-submitting to avoid redundant computation.

## Example: VADR Batch Analysis

From the GATK-to-FreeBayes regression testing (PR #1053), we ran VADR on 30 FASTAs
(15 assemblies x old/new) using dsub:

```bash
source ~/venvs/dsub/bin/activate  # or wherever dsub is installed

dsub --provider google-cls-v2 \
  --project gcid-viral-seq \
  --regions us-central1 \
  --image staphb/vadr:1.6.4 \
  --machine-type n1-highmem-4 \
  --script run_vadr.sh \
  --tasks vadr_tasks.tsv \
  --logging gs://viral-temp-30d/vadr_regression/logs/
```

The tasks TSV had columns for VADR options (`--env VADR_OPTS`), model URL
(`--env MODEL_URL`), input FASTA (`--input FASTA`), and outputs
(`--output NUM_ALERTS`, `--output ALERTS_TSV`, `--output VADR_TGZ`).

All 30 jobs completed in ~15 minutes total (running in parallel on GCP).
