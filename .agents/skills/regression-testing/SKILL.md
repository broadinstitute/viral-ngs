# Assembly Regression Testing

End-to-end regression testing for assembly pipeline changes against Terra submissions.

## When to Use

Use this playbook when a PR makes functional changes to the assembly or variant-calling
pipeline (e.g., swapping variant callers, changing alignment parameters, modifying
consensus logic). It compares assembly outputs from old vs new code across hundreds
of real samples to validate equivalence or improvement.

## Prerequisites

- **gcloud CLI** -- authenticated with access to Terra workspace GCS buckets
- **mafft** -- for pairwise sequence alignment
- **Python** with pandas and matplotlib (e.g., a dataviz venv)
- **dsub** -- for running VADR batch jobs on GCP (see the `dsub-batch-jobs` skill)

## Workflow

### Step 1: Set Up Terra Submissions (Manual)

The user must manually launch Terra submissions with old and new code:
1. Run the pipeline on a representative dataset using the **main branch** Docker image
2. Run the same pipeline on the same dataset using the **feature branch** Docker image
3. Note the submission IDs and workspace bucket for both runs

### Step 2: Discover Paired Samples

Use `discover_pairs.py` to find all comparable old/new sample pairs by crawling
GCS Cromwell output directories.

```bash
python discover_pairs.py <workspace_name> \
  --bucket <workspace-bucket-id> \
  --old-sub <old-submission-id> \
  --new-sub <new-submission-id> \
  -o pairs.json
```

This produces a JSON mapping sample_name -> {old_tsv, new_tsv} for all samples
present in both submissions.

### Step 3: Compare Assembly Outputs

Use `compare_sample_pair.py` to compare each sample pair. This script:
- Downloads assembly_stats TSVs from GCS
- Compares metrics (coverage, % reference covered, length, etc.)
- Downloads FASTA assemblies and aligns them with mafft
- Reports SNPs, indels (events and bp), ambiguity diffs, and terminal extensions

```bash
python compare_sample_pair.py \
  --old-tsv <gcs_uri> --new-tsv <gcs_uri> \
  --work-dir ./results/<sample> \
  --output-json ./results/<sample>.json
```

For batch processing, use `run_regression.py` (in the terra-regression/scripts/ directory)
which orchestrates parallel execution across all pairs.

### Step 4: Generate Report

Use `generate_report.py` to aggregate all per-sample JSONs into plots and a markdown report.

```bash
python generate_report.py \
  --results-dir ./results/ \
  --report-dir ./report/ \
  --workspace-name <name>
```

Outputs:
- Summary TSV with per-assembly metrics
- 8 plots (scatter plots, histograms, identity distribution)
- Markdown report with summary tables and divergent assembly details

### Step 5: (Optional) VADR Annotation Quality

For assemblies with internal indel differences, run VADR to assess whether indels
cause frameshifts or other annotation problems. See the `dsub-batch-jobs` skill for
details on running batch jobs via dsub.

Use `run_vadr.sh` with dsub to run VADR on each FASTA:

```bash
dsub --provider google-cls-v2 \
  --project <gcp-project> --regions us-central1 \
  --image staphb/vadr:1.6.4 \
  --machine-type n1-highmem-4 \
  --script run_vadr.sh \
  --tasks vadr_tasks.tsv \
  --logging gs://<bucket>/vadr_logs/
```

VADR model parameters come from the viral-references repo:
https://github.com/broadinstitute/viral-references/blob/main/annotation/vadr/vadr-by-taxid.tsv

Use the taxid from the assembly_id (format: `sample_id-taxid`) to look up the
correct `vadr_opts`, `min_seq_len`, `max_seq_len`, `vadr_model_tar_url`, and
`vadr_model_tar_subdir`.

### Step 6: Post Results

Post the report as a PR comment. Before posting:
- **Self-review the proposed comment for confidential information** (sample names,
  internal paths, credentials, etc.). Ask the user if in doubt about what is safe
  to share publicly.
- Include plots as image attachments if the PR is on GitHub
- Attribute the analysis appropriately

## Key Patterns

### Per-Segment Alignment

Multi-segment genomes (e.g., influenza with 8 segments) must be aligned
**per-segment**, not as a single concatenated sequence. Otherwise, terminal
effects at segment boundaries get misclassified as internal indels.

The `compare_sample_pair.py` script handles this automatically: it pairs
segments by FASTA header, aligns each pair independently, analyzes each
alignment separately (so terminal effects stay terminal), and aggregates
the statistics.

### Event Counting vs BP Counting

For indels, both counts matter:
- **BP count**: Total gap positions (e.g., "49 bp of indels")
- **Event count**: Contiguous gap runs (e.g., "13 indel events")

A single 26bp insertion is 1 event but 26 bp. Event counts better reflect
the number of variant-calling decisions that differ between old and new code.

### VADR Frameshift Cascade Detection

A single spurious 1bp indel in a coding region causes a cascade of VADR alerts:
1. `FRAMESHIFT` -- the indel shifts the reading frame
2. `STOP_CODON` -- premature stop codon in the shifted frame
3. `UNEXPECTED_LENGTH` -- protein length doesn't match model
4. `PEPTIDE_TRANSLATION_PROBLEM` -- for each downstream mature peptide

When comparing VADR alert counts, a large delta (e.g., 32 -> 1) usually means
one version has frameshift-causing indels that the other avoids. Check the
`.alt.list` files to confirm which genes are affected.

## Interpreting Results

### What to Look For

1. **Identity distribution**: Most assemblies should be 100% identical. Any
   below 99.9% warrant investigation.
2. **SNP count = 0 for all assemblies**: Pipeline changes that only affect
   indel calling (e.g., swapping variant callers) should produce zero SNPs.
3. **Indel events**: The number and nature of indel differences. Are they in
   coding regions? Do they cause frameshifts?
4. **Coverage correlation**: Low-coverage samples (<10x) are most likely to
   show differences between variant callers.
5. **VADR alert deltas**: Fewer alerts = more biologically plausible assembly.
   Large improvements (e.g., -31 alerts) strongly favor the new code.

### Red Flags

- Assemblies present in old but missing in new (or vice versa)
- SNPs introduced where none existed before
- VADR alerts increasing significantly for the new code
- Differences concentrated in specific organisms/taxids
