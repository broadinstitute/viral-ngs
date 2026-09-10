# VirNucPro post-processing

This document describes the post-processing logic that `viral_ngs.classify.virnucpro` layers on top of the raw VirNucPro classifier output. It is purely descriptive: it explains what the code does today so that reviewers can evaluate the choices independently. It does **not** prescribe what the logic should be.

> **Scope.** The VirNucPro classifier itself (a GPU/LLM-based viral-vs-non-viral predictor) lives in a separate Docker image and is not part of this repository. This module only handles CPU-friendly post-processing of VirNucPro outputs.

---

## Pipeline overview

```
    contigs.fasta                                reads.fastq
         │                                            │
         ▼                                            ▼
   ┌─────────────┐                            ┌──────────────┐
   │  VirNucPro  │  (separate Docker image)   │   Aligner    │  (any aligner; minimap2, BWA, …
   │ (GPU/LLM)   │                            │              │   anything that emits NM tag)
   └─────────────┘                            └──────────────┘
         │                                            │
         │  highest_scores.tsv                        │  aligned.bam
         │  (one row per chunk)                       │  (reads aligned to contigs)
         ▼                                            │
   ┌─────────────────────────────┐                    │
   │   classify_contigs()        │  ◄────── Stage 1   │
   │   (this module)             │                    │
   └─────────────────────────────┘                    │
         │                                            │
         │  contigs.tsv                               │
         │  (one row per contig:                      │
         │   call + tier + metrics)                   │
         │                                            │
         └─────────────┬──────────────────────────────┘
                       ▼
              ┌──────────────────────────────────┐
              │   classify_reads_by_contig()     │  ◄────── Stage 2
              │   (this module)                  │
              └──────────────────────────────────┘
                       │
                       ▼
                 reads_classified.tsv
                 (one row per query_name from the BAM)
```

Two stages, two CLI subcommands:

| CLI subcommand                       | Helper                          | Stage                          |
|--------------------------------------|---------------------------------|--------------------------------|
| `virnucpro_contigs`                  | `classify_contigs`              | Chunk → contig roll-up         |
| `virnucpro_label_reads_by_contig`    | `classify_reads_by_contig`      | Contig → read label propagation|

---

## Background: what VirNucPro produces

VirNucPro chops each input contig into fixed-length non-overlapping chunks (the published benchmarking work uses 300 bp and 500 bp), runs a six-frame translation, scores each frame with an LLM, and emits two TSVs. We consume `highest_scores.tsv`, which has one row per chunk:

| column         | meaning                                                       |
|----------------|---------------------------------------------------------------|
| `Modified_ID`  | `{contig_id}_chunk_{N}` — N is sparse (chunks failing internal VirNucPro filters are omitted) |
| `Is_Virus`     | boolean; empirically equal to `max_score_1 > max_score_0`     |
| `max_score_0`  | non-viral class score (probability-like, 0–1)                 |
| `max_score_1`  | viral class score (probability-like, 0–1)                     |

`Is_Virus` is not consumed by this module — the two score columns carry strictly more information.

---

## Stage 1: chunk → contig (`classify_contigs`)

Takes a VirNucPro `highest_scores.tsv` and produces a per-contig TSV with a categorical `call`, a `tier` confidence label, and supporting metrics.

### Grouping chunks back into contigs

VirNucPro flattens the contig structure into chunk IDs like `NODE_1_length_15270_cov_X_chunk_5`. To re-aggregate, the helper does two things to the `Modified_ID` column:

1. **Contig label** — strip the trailing `_chunk_\d+` suffix. Used as the human-readable `ID` in the output.
2. **Group key** — match against `id_pattern` (default `(NODE_\d+)`) and use the matched group as the grouping key for `pandas.DataFrame.groupby`.

These two derivations are independent. With the default SPAdes-style pattern, the group key is the `NODE_N` prefix and the contig label is the full pre-chunk string (`NODE_N_length_X_cov_Y`). Rows whose ID does not match `id_pattern` are dropped with a warning. If **no** rows match, `ValueError` is raised. If the input file is empty or the table is header-only, an empty (header-only) output TSV is written and the helper returns without error.

### Per-chunk derived quantities

For each chunk in a contig group:

```
delta      = max_score_1 - max_score_0          # signed; positive = viral, negative = non-viral
confidence = sqrt(|delta|)                       # confidence weight for the per-contig average
```

Three orthogonal categorical labels are also computed at the chunk level, used later for tier assignment:

```
confident_viral    = (max_score_1 > min_confident_score) AND (max_score_0 < max_opposing_score)
confident_nonviral = (max_score_0 > min_confident_score) AND (max_score_1 < max_opposing_score)
ambiguous          = (max_score_1 > min_ambiguous_score) AND (max_score_0 > min_ambiguous_score)
```

A chunk can be **none** of the three (e.g. `max_score_1 = 0.6, max_score_0 = 0.4`), exactly one, or `ambiguous` plus possibly one of the confident flags (depending on parameter values). The default parameters (`min_confident_score=0.8`, `max_opposing_score=0.3`, `min_ambiguous_score=0.7`) make these flags effectively mutually exclusive.

### Per-contig aggregated quantities

```
weighted_delta = Σ(delta_i * confidence_i) / Σ confidence_i
                 (falls back to mean(delta) if all confidences are zero)

n_chunks                = total chunks in the group
n_confident_viral       = count of chunks where confident_viral is true
n_confident_nonviral    = count of chunks where confident_nonviral is true
n_ambiguous             = count of chunks where ambiguous is true
n_effective             = n_chunks - n_ambiguous
viral_proportion        = n_confident_viral    / n_effective    (0 if n_effective == 0)
nonviral_proportion     = n_confident_nonviral / n_effective    (0 if n_effective == 0)
```

### Call decision

A single threshold splits the three call categories:

```
if   weighted_delta >  min_weighted_delta : call = "Viral"
elif weighted_delta < -min_weighted_delta : call = "Non-viral"
else                                       : call = "Ambiguous"
```

So `call` is fully determined by `weighted_delta` and `min_weighted_delta` (default 0.3). The other parameters do not move contigs between call categories — they only affect `tier`.

### Tier decision

```
Ambiguous call:
    tier = "review"

Viral call:
    if n_confident_viral >= 1 AND viral_proportion >= min_viral_proportion:
        tier = "high_confidence"  if weighted_delta >  high_confidence_delta else "moderate_confidence"
    else:
        tier = "low_confidence"

Non-viral call:  (mirror of Viral; uses nonviral_proportion and the negative delta)
    if n_confident_nonviral >= 1 AND nonviral_proportion >= min_nonviral_proportion:
        tier = "high_confidence"  if weighted_delta < -high_confidence_delta else "moderate_confidence"
    else:
        tier = "low_confidence"
```

Finally, a **short-contig demotion** is applied:

```
if n_chunks < min_chunks:
    "high_confidence" or "moderate_confidence"  →  "low_confidence"
    "low_confidence"                            →  "review"
    # "review" is unchanged
```

### Parameters

| Parameter (helper kwarg / CLI flag)          | Default | Affects | Meaning |
|---------------------------------------------|---------|---------|---------|
| `min_weighted_delta` / `--min-weighted-delta` | 0.3     | **call** + tier | Min `|weighted_delta|` for a Viral/Non-viral call. Below: Ambiguous. |
| `high_confidence_delta` / `--high-confidence-delta` | 0.6 | tier | Split between `moderate_confidence` and `high_confidence`. |
| `min_confident_score` / `--min-confident-score` | 0.8  | tier | Winning-class score for `confident_viral`/`confident_nonviral` flags. |
| `max_opposing_score` / `--max-opposing-score` | 0.3   | tier | Max opposing-class score for the same. |
| `min_ambiguous_score` / `--min-ambiguous-score` | 0.7 | tier | Min score in **both** classes for the `ambiguous` flag (and `n_effective` denominator). |
| `min_viral_prop` / `--min-viral-prop`        | 0.1     | tier | Min `viral_proportion` for non-`low_confidence` Viral tiers. |
| `min_nonviral_prop` / `--min-nonviral-prop`  | 0.1     | tier | Mirror for Non-viral. |
| `min_chunks` / `--min-chunks`                | 5       | tier | Below this, tier is demoted one notch. |
| `id_col` / `--id-col`                        | `"Modified_ID"` | grouping | Source column for chunk IDs. |
| `id_pattern` / `--id-pattern`                | `(NODE_\d+)` | grouping | Regex group used to roll chunks back into contigs. |

### Output schema (`classify_contigs`)

| column                 | type   | notes |
|------------------------|--------|-------|
| `ID`                   | string | chunk-suffix-stripped contig label (e.g. `NODE_1_length_15270_cov_X`) |
| `call`                 | string | one of `Viral`, `Non-viral`, `Ambiguous` |
| `tier`                 | string | one of `high_confidence`, `moderate_confidence`, `low_confidence`, `review` |
| `weighted_delta`       | float  | rounded to 3 decimals |
| `n_chunks`             | int    | |
| `n_confident_viral`    | int    | |
| `n_confident_nonviral` | int    | |
| `n_ambiguous`          | int    | |
| `viral_proportion`     | float  | rounded to 3 decimals |
| `nonviral_proportion`  | float  | rounded to 3 decimals |

Rows are sorted by `ID` with a natural-sort key (so `NODE_2` precedes `NODE_10`).

---

## Stage 2: contig → read (`classify_reads_by_contig`)

Takes (a) an aligned BAM of reads-to-contigs, and (b) the per-contig TSV from Stage 1, and emits a per-read TSV with each read's call/tier propagated from its best-mapping contig.

### Step A: BAM filtering and per-alignment augmentation

`_prepare_augmented_bam_file` iterates the BAM with `pysam.AlignmentFile` and emits one row per **retained** alignment to a normalized intermediate TSV. Rows are filtered out **before** they reach DuckDB:

- `record.is_unmapped` (flag `0x4`)              → dropped
- `record.is_secondary` (flag `0x100`)           → dropped (counted as `n_secondary`)
- `record.is_supplementary` (flag `0x800`)       → dropped (counted as `n_secondary`)

Everything else (i.e. **primary mapped alignments**) is emitted. The flow does not inspect the paired-end flags (`0x1`, `0x40`, `0x80`); R1 and R2 of a pair both pass the filter if both are primary mapped, and they appear as two separate rows that share `query_name`.

For each retained alignment, the helper validates and computes:

| field                  | derivation |
|------------------------|------------|
| `source_order`         | 1-based counter assigned in BAM iteration order (final tiebreaker for best-alignment selection) |
| `query_name`           | `record.query_name` |
| `query_length`         | `record.infer_read_length()` → `infer_query_length(always=True)` → `record.query_length` (first non-None) |
| `query_start`/`_end`   | `record.query_alignment_start` / `record.query_alignment_end` (excludes soft clips) |
| `strand`               | `'-'` if `record.is_reverse` else `'+'` |
| `target_name`          | reference name from `bam.get_reference_name(record.reference_id)` |
| `target_length`        | `bam.get_reference_length(target_name)` |
| `target_start`/`_end`  | `record.reference_start` / `record.reference_end` |
| `alignment_block_length` | sum of CIGAR op lengths for `M`, `I`, `D`, `=`, `X` |
| `num_matches`          | `alignment_block_length − NM` |
| `pct_identity`         | `100 × num_matches / alignment_block_length` (BLAST-like; indels in denominator) |
| `pct_query_cov`        | `100 × (query_end − query_start) / query_length` |
| `mapping_quality`      | `record.mapping_quality` |

Validation errors that raise `ValueError`:

- Missing `NM` tag on any retained primary alignment.
- `query_length` cannot be inferred or is zero.
- `reference_end` cannot be inferred.
- `alignment_block_length == 0`.
- `NM > alignment_block_length`.

### Step B: DuckDB pipeline

The normalized alignment TSV and a normalized copy of the contig classifications are loaded into an in-memory DuckDB connection. DuckDB memory is capped from the cgroup limit (`~75%` of the container's memory) unless overridden by `--duckdb-memory-limit`; passing an empty string opts out of any limit.

A `mapped_well` flag is computed but **not used to filter rows**:

```sql
mapped_well = (mapping_quality >= min_mapq
               AND pct_identity   >= min_identity
               AND pct_query_cov  >= min_query_cov)
```

The retained-primary alignments are LEFT-JOINed to the contig classifications. Contigs absent from the classifications table get a synthesized call:

```sql
COALESCE(c.call, 'Unclassified')   AS call
COALESCE(c.tier, '')               AS tier
COALESCE(c.<other_metrics>, 0)     AS <other_metrics>
```

This `merged` table has one row per retained primary alignment, with the contig's call attached (or `Unclassified`).

Two derived tables are then computed in parallel from `merged`:

**`call_counts`** — number of distinct call categories per `query_name`, across **all** rows in `merged`:

```sql
SELECT query_name, COUNT(DISTINCT call) AS n_distinct_calls
FROM merged
GROUP BY query_name
```

Notes on what this counts:
- It groups by `query_name` only. R1 and R2 of a paired-end pair share `query_name`, so their alignments both contribute.
- It treats `'Unclassified'` as a value alongside `'Viral'`, `'Non-viral'`, `'Ambiguous'` — so a read whose alignments touch one classified contig and one unclassified contig has `n_distinct_calls == 2`.
- It does **not** filter by `mapped_well` — sub-threshold alignments contribute to the count just like above-threshold ones.

**`best`** — single winning alignment per `query_name`:

```sql
ROW_NUMBER() OVER (
    PARTITION BY query_name
    ORDER BY mapping_quality DESC, pct_identity DESC, source_order ASC
) AS rn
WHERE rn = 1
```

Notes:
- Partitioning by `query_name` collapses paired-end mates (and any other case where one `query_name` has multiple primary alignments) down to a single row.
- The losing alignment(s) are dropped from the output entirely — they are not represented as separate rows.
- `source_order` (BAM iteration index) is the final tiebreaker, so two alignments with identical `mapping_quality` and `pct_identity` are resolved by which appeared first in the input BAM. Reordering the BAM between alignment and classification can change which mate "wins" in true ties.

### Step C: result composition

The final result joins `best` to `call_counts` and applies a `Multi-mapped` override:

```sql
call = CASE WHEN cc.n_distinct_calls > 1 THEN 'Multi-mapped' ELSE b.call END
tier = CASE WHEN cc.n_distinct_calls > 1 THEN 'review'       ELSE b.tier END
```

So the per-read `call` can take five values:

| value         | source |
|---------------|--------|
| `Viral`       | propagated from the winning contig's Stage-1 call |
| `Non-viral`   | propagated from the winning contig's Stage-1 call |
| `Ambiguous`   | propagated from the winning contig's Stage-1 call |
| `Unclassified`| synthesized when the winning contig is missing from the classifications table |
| `Multi-mapped`| synthesized when `n_distinct_calls > 1` for the read's `query_name` (overrides the four above) |

### Output schema (`classify_reads_by_contig`)

One row per `query_name` that had at least one retained primary alignment. Unmapped reads, secondary alignments, and supplementary alignments produce no output rows.

| column                 | source |
|------------------------|--------|
| `read_id`              | `query_name` of the winning alignment |
| `read_length`          | `query_length` of the winning alignment |
| `contig_id`            | `target_name` of the winning alignment |
| `contig_length`        | `target_length` of the winning alignment |
| `strand`               | from the winning alignment |
| `mapping_quality`      | from the winning alignment |
| `pct_identity`         | from the winning alignment |
| `pct_query_cov`        | from the winning alignment |
| `mapped_well`          | from the winning alignment (boolean) |
| `call`                 | one of the five values above |
| `tier`                 | from Stage 1, or `'review'` when overridden by Multi-mapped, or `''` for Unclassified |
| `weighted_delta`       | from Stage 1 (0 for Unclassified) |
| `n_chunks`             | from Stage 1 (0 for Unclassified) |
| `n_confident_viral`    | from Stage 1 (0 for Unclassified) |
| `n_confident_nonviral` | from Stage 1 (0 for Unclassified) |
| `n_ambiguous`          | from Stage 1 (0 for Unclassified) |
| `viral_proportion`     | from Stage 1 (0 for Unclassified) |
| `nonviral_proportion`  | from Stage 1 (0 for Unclassified) |

Output row order is whatever DuckDB's `COPY` emits (no explicit `ORDER BY`). The output extension governs compression: `.gz`, `.zst`, `.lz4`, `.bz2` are detected automatically; anything else is written uncompressed.

If the BAM yields zero retained primary alignments (empty BAM or all reads unmapped/secondary/supplementary), an empty header-only output TSV is written and the DuckDB pipeline is skipped.

### Parameters

| Parameter (CLI flag)                | Default     | Meaning |
|------------------------------------|-------------|---------|
| `--min-mapq`                       | 5           | Threshold for the `mapped_well` flag. Sub-threshold alignments are **not** dropped. |
| `--min-identity`                   | 90.0        | Percent identity threshold for `mapped_well`. Percent units (90, not 0.9); a fractional value 0 < x < 1 raises `ValueError`. |
| `--min-query-cov`                  | 80.0        | Percent query coverage threshold for `mapped_well`. Same percent-units rule. |
| `--duckdb-memory-limit`            | auto-detect | DuckDB memory cap (e.g. `"8GB"`). Default auto-detects ~75% of the cgroup limit. Empty string opts out of any cap. |
| `--work-dir`                       | system tmp  | Parent directory for the per-run scratch tree (used for the normalized TSVs and DuckDB spill). Distinct from the `--tmp_dir` flag provided by `cmd.common_args`. |

---

## Behavioral notes for reviewers

The following are factual statements about the implementation, not recommendations.

1. **Primary-alignment filter is at the pysam stage, not in SQL.** Secondary (`0x100`) and supplementary (`0x800`) alignments never reach DuckDB. The `ROW_NUMBER()` window therefore exists only to deduplicate `query_name` values whose primary alignments appear more than once in the BAM — i.e. paired-end mates or pathological multi-primary cases.

2. **Paired-end mates share `query_name` and are not distinguished.** Both the `PARTITION BY query_name` (winner selection) and the `GROUP BY query_name` (distinct-call count) treat R1 and R2 of a pair as alignments of the same logical read. There is no use of the BAM paired-end flags (`0x1`, `0x40`, `0x80`).

3. **`Unclassified` participates in the `Multi-mapped` count.** Because `COALESCE(c.call, 'Unclassified')` is applied before `COUNT(DISTINCT call)`, a `query_name` with one alignment to a classified contig and one to an unclassified contig has `n_distinct_calls == 2` and is labeled `Multi-mapped`.

4. **`mapped_well` is informational only.** It is computed and exposed as an output column, but no DuckDB query filters on it. Sub-threshold alignments contribute to both `best` selection (where high-quality alignments will still win on `mapping_quality DESC` ordering) and `call_counts` (where any call value they pull from the classifications table is counted as a distinct call).

5. **Output row ordering is not guaranteed.** No explicit `ORDER BY` on the terminal `COPY`.

6. **Default thresholds in percent units.** `min_identity` and `min_query_cov` validate against `0 < value < 1.0` and raise — so fractional thresholds cannot be silently passed. Zero and ≥1.0 are accepted.
