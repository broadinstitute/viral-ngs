# Kallisto Implementation — Code Review Findings

Source review conducted while writing the Kallisto section of
`src/viral_ngs/classify/DESIGN.md`. All findings are **pre-existing** in the
current branch and are independent of the VirNucPro CLI work and the DESIGN.md
documentation edits.

Files in scope:
- `src/viral_ngs/classify/kallisto.py`
- `src/viral_ngs/metagenomics.py` (kallisto CLI entry points)

Status legend: **[correctness]** produces silent wrong/missing output in
production · **[maintainability]** code health · **[verify]** suspected issue,
not yet confirmed.

---

## 1. Missing input file silently succeeds — [correctness]

`kallisto.py:194-196` (`classify`) and `kallisto.py:277-278` (`extract`):

```python
if self._is_fastq(in_bam):
    if not os.path.exists(in_bam):
        return
```

A non-existent FASTQ input makes the command exit 0 with **no outputs and no
error**. This hides pipeline misconfiguration: a missing input should raise a
loud `FileNotFoundError`, not be treated as a successful no-op.

This is the behavior behind the FASTQ caveat now documented in DESIGN.md
("Empty Input Contract"). The doc describes the behavior; the behavior itself
is the defect.

**Recommendation:** raise `FileNotFoundError(in_bam)` when the input path does
not exist, in both `classify` and `extract`.

**Proposed solution:**

- Replace both FASTQ missing-file `return` branches with
  `raise FileNotFoundError(in_bam)`.
- Add focused unit tests for `classify()` and `extract()` with a missing FASTQ
  path, asserting that `FileNotFoundError` is raised.
- Update `src/viral_ngs/classify/DESIGN.md` after the code change to remove the
  current caveat that a missing FASTQ silently returns without outputs.

---

## 2. Empty extract input writes nothing — asymmetric with count — [correctness]

`kallisto.py:283-284` (`extract`):

```python
else:
    if samtools.SamtoolsTool().isEmpty(in_bam):
        return
```

Count handles empty BAM by writing stable empty artifacts via
`_write_empty_count_outputs` (`kallisto.py:359-369`): a `1 x 0` `adata.h5ad`
plus a header-only `counts.tsv`. Extract, by contrast, returns early and writes
**neither `summary.tsv` nor `read_hits.tsv`**.

Any WDL/pipeline task that declares `summary.tsv` (or `read_hits.tsv`) as a
required output will break on an empty sample, even though the "stable
artifacts for empty input" philosophy is explicitly part of the count contract.

**Recommendation:** add a parallel `_write_empty_extract_outputs` that emits
header-only `summary.tsv` and `read_hits.tsv`, mirroring the count contract, and
call it from the empty-BAM branch of `extract`.

**Proposed solution:**

- Add `_write_empty_extract_outputs(out_dir)` to create `out_dir` and write:
  - `read_hits.tsv` with header `read_id	db_hit_id`
  - `summary.tsv` with header `SAMPLE_ID	READ_ID	DB_ID	TAXONOMY_LINEAGE	TAXONOMY_NAME	SEQUENCE_LENGTH`
- In `extract()`, call `_write_empty_extract_outputs(out_dir)` before returning
  from the empty-BAM branch.
- Add a behavior-level unit test that runs `extract()` with `SamtoolsTool.isEmpty`
  mocked true and asserts the two header-only TSVs exist. Avoid asserting private
  helper call counts as the primary test behavior.
- Add a direct helper test only if the header construction is not already fully
  covered by the behavior-level test.

---

## 3. Multi-sample H5AD silently sums across samples — [correctness]

`extract_hit_ids_from_h5ad` (`kallisto.py:712-732`) relies on `_h5ad_hit_totals`
(`kallisto.py:734-741`):

```python
gene_totals = adata.X.sum(axis=0)   # sums across ALL observations/samples
```

The documented and intended contract is single-sample (see the method docstring
"Assumes h5ad contains a single sample" and the acknowledging TODO at
`metagenomics.py:1428`). A multi-sample H5AD does not error — it silently sums
counts across samples and returns the union of targets, which is almost
certainly not what the caller wants.

**Recommendation:** assert `adata.n_obs == 1` and raise a clear error otherwise,
matching the documented assumption. This removes a silent-wrong-result path and
lets the existing TODO at `metagenomics.py:1428` be deleted.

**Proposed solution:**

- In `extract_hit_ids_from_h5ad()`, after `anndata.read_h5ad()`, validate
  `adata.n_obs == 1`.
- Raise a clear `ValueError` such as
  `Expected single-sample h5ad for kallisto extract target selection; found {adata.n_obs} observations`.
- Add a unit test with a two-row H5AD showing that target extraction raises
  instead of summing across samples.
- Remove the TODO in `metagenomics.py` once the assumption is enforced in
  `extract_hit_ids_from_h5ad()`.

---

## 4. `execute()` detects failure by grepping stderr text — [maintainability]

`kallisto.py:121-131`:

```python
has_error = stderr and ('Traceback (most recent call last)' in stderr or 'CalledProcessError' in stderr)
if process.returncode != 0 or has_error:
    ...
    raise subprocess.CalledProcessError(...)
```

The comment notes "the kb backend sometimes catches exceptions without proper
exit codes," so return code alone is unreliable. The substring scan is fragile
in both directions:
- **false positive:** incidental occurrence of those strings in stderr
- **false negative:** kb swallowing an error in a form that prints neither
  marker

Fully fixing this needs kb-python cooperation / version pinning, so it is not a
quick change.

**Recommendation:** track as a follow-up. At minimum, add a focused test
reproducing the "nonzero work, zero exit code" case this guards against, and a
comment naming the exact kb-python version(s) and failure mode it compensates
for, so it can be removed when kb behavior is fixed.

**Proposed solution:**

- Do not remove the stderr guard in this PR unless we can prove current
  `kb-python` returns nonzero for the historical failure mode.
- Add a small unit test for `execute()` where mocked `Popen` returns code `0`
  but stderr contains a traceback marker, asserting that
  `subprocess.CalledProcessError` is raised. This documents the intentional
  behavior.
- Add a TODO/comment with the observed `kb-python` version and command shape
  that motivated the guard if we can reproduce it during testing.
- Track a follow-up to remove the stderr substring detection once the wrapper can
  rely on `kb-python` exit codes.

---

## 5. Dead / test-only methods — [maintainability]

- `parse_h5ad_counts` (`kallisto.py:700-710`): **no callers anywhere** in `src`
  or `tests`. Delete.
- `_write_read_hits_tsv` (`kallisto.py:502-503`): referenced **only** by tests
  (`tests/unit/classify/test_tools_kallisto.py:449,498`), never by production
  code.
- `_iter_fastq_read_ids` (`kallisto.py:650-653`): referenced **only** by tests
  (`test_tools_kallisto.py:506,512`), never by production code.

Testing private helpers that nothing calls is coverage theater.

**Recommendation:** delete `parse_h5ad_counts`. For `_write_read_hits_tsv` and
`_iter_fastq_read_ids`, either wire them into the real path or remove them and
fold their assertions into tests of the actual entry points (`extract` /
`_write_extract_tsvs` / `_iter_fastq_records`).

**Proposed solution:**

- Delete `parse_h5ad_counts()`.
- Delete `_write_read_hits_tsv()` and update tests currently using it to call
  `_write_extract_tsvs()` instead, because `_write_extract_tsvs()` is the real
  production path.
- Delete `_iter_fastq_read_ids()` and update malformed FASTQ tests to exercise
  `_iter_fastq_records()` directly.
- Keep tests focused on observable output files and parser validation rather
  than private wrappers that are not called by production code.

---

## 6. `execute()` docstring drift — [maintainability]

`kallisto.py:78-86`: the docstring documents an `input:` parameter that is not
in the signature `execute(self, command, output, args=None, options=None)`.

**Recommendation:** remove the stale `input:` line from the docstring.

**Proposed solution:**

- Remove the `input:` entry from the `execute()` docstring.
- While touching the docstring, keep the existing `command`, `output`, `args`,
  and `options` descriptions aligned with the actual signature.

---

## 7. Single-cell `technology` choices vs. generic BAM interleaving — [verify]

`parser_kallisto` exposes single-cell technologies
(`metagenomics.py:397-399`: `10xv2`, `10xv3`, `dropseq`, ...), but the BAM→FASTQ
ingestion path runs Picard `SamToFastq` plus naive R1/R2 interleaving
(`kallisto.py:220-236`). It is unverified whether cell barcodes survive that
path.

**Not asserted as a bug** — needs a barcode round-trip test to confirm. If only
`bulk` is actually supported end-to-end from BAM input, the single-cell choices
list is a trap and should be removed or gated.

**Recommendation:** add a barcode-preservation test for one single-cell
technology from BAM input; if it does not survive, restrict the `--technology`
choices accordingly.

**Proposed solution:**

- Treat as a tracked follow-up, not a blocker for the TSV/refactor PR, unless the
  current pipeline intends to run single-cell Kallisto modes from BAM immediately.
- Add an integration-style fixture with a small BAM carrying the barcode/UMI tags
  needed by one supported single-cell technology.
- Run the BAM through the current Picard `SamToFastq` + interleaving path and
  verify whether `kb count -x <technology>` receives enough information to
  preserve barcode semantics.
- If the test fails, restrict `--technology` for BAM inputs to `bulk` or require
  FASTQ input for single-cell technologies, and document that constraint in
  `DESIGN.md`.

---

## Suggested batching

| Group | Findings | Notes |
| --- | --- | --- |
| Correctness PR (with tests) | 1, 2, 3 | Each produces silent wrong/missing output in production. |
| Trivial cleanup (same files) | 5, 6 | Fold into the correctness PR. |
| Tracked follow-ups | 4, 7 | 4 needs kb-version investigation; 7 needs a barcode round-trip test. |

---

# Review of Proposed Solutions

Feedback on the **Proposed solution** blocks above, after re-verifying the
load-bearing assumptions against the source. Conclusions are grounded in the
code as it currently exists on this branch.

## Blocking issue: Findings 1 + 2 ship an incoherent contract as proposed

`samtools.SamtoolsTool().isEmpty` (`samtools.py:469-471`):

```python
def isEmpty(self, inBam):
    if not os.path.isfile(inBam):
        return True
```

**`isEmpty` returns `True` for a file that does not exist.** A missing *BAM*
therefore never reaches the FASTQ existence check in Finding 1 — it flows into
the `else` branch, `isEmpty` reports `True`, and the code takes the empty-input
path. Applying both proposals as written yields:

| Input | Missing-file behavior after proposed fixes |
| --- | --- |
| FASTQ | `FileNotFoundError` (loud) — correct |
| BAM | silently treated as "empty", writes empty outputs — WRONG |

BAM is the primary documented input (`parser.add_argument('in_bam', ...)`), so
the proposal raises on the rare case and leaves the common case silently
succeeding — the exact failure mode Finding 1 claims to eliminate. Finding 1's
writeup only cites the FASTQ branch; neither the finding nor its proposed
solution noticed that count/extract already swallow a missing BAM via `isEmpty`.

Worse, with Finding 2 applied, a typo'd BAM path would now *write empty extract
outputs* — upgrading a silent no-op into silent fake-success.

**Required rework:** treat Findings 1 and 2 as one input-handling contract, not
two independent patches. Hoist a single existence check above the FASTQ/BAM
branch in both `classify` and `extract`:

```python
if not os.path.exists(in_bam):
    raise FileNotFoundError(in_bam)
```

Only then does the `else` branch's `isEmpty` mean "genuinely empty BAM," making
Finding 2's empty-output path reachable only for real empty input. Confirmed no
existing test passes a missing FASTQ expecting a silent return
(`test_tools_kallisto.py` only mocks `isEmpty`), so this does not break current
tests.

## Per-finding verdicts

### Finding 1 — direction right, scope too narrow
- Must cover BAM, not just FASTQ (see blocking issue). Hoist the check above the
  branch.
- DESIGN.md follow-up should *rewrite* the empty-input section to state "missing
  input (BAM or FASTQ) raises," not merely delete the FASTQ caveat.

### Finding 2 — sound, two refinements
- **DRY the headers.** The proposed `_write_empty_extract_outputs` re-hardcodes
  the `read_hits.tsv` / `summary.tsv` headers that `_write_extract_tsvs` already
  writes (`kallisto.py:515-516`); a future column change would drift across two
  sites. Share a header constant/tuple between both paths (keep the existing
  empty-`target_ids` guard at `kallisto.py:506-507`).
- **Unaddressed asymmetry.** In `--h5ad` derive-targets mode, an empty-count
  H5AD makes `extract_hit_ids_from_h5ad` return `[]`, and `kallisto_extract`
  then raises `ValueError` (`metagenomics.py:1432-1433`) before `extract()`
  runs — so extract still hard-errors on an empty sample in that mode and never
  reaches the new empty-output path. The proposal only stabilizes the explicit
  `--targets` + empty-BAM path. Decide deliberately whether empty-via-derive
  should also produce stable empty artifacts or stay an error, and record the
  decision in DESIGN.md (the proposal does not mention a doc update for #2).

### Finding 3 — correct, no conflicts
- Verified the empty-count H5AD is `1 x 0` so `n_obs == 1`: the new assert
  passes and still returns `[]` cleanly. Single caller is
  `metagenomics.py:1430`. Ready to implement as written.

### Finding 4 — good, conservative
- Keeping the guard plus a regression test ("exit 0 + traceback in stderr →
  raises") is the right call. The test must patch `subprocess.Popen` in the
  `kallisto` module namespace and set `communicate.return_value` +
  `returncode = 0`; the existing mocking pattern at `test_tools_kallisto.py:120-121`
  can be reused.

### Finding 5 — sound; two notes
- `parse_h5ad_counts` is nominally **public** (no leading underscore), so its
  deletion is technically an API change. Zero in-repo callers confirmed and it
  is undocumented, so risk is low — call it out in the PR description in case
  viral-pipelines reaches into it.
- Test migrations are behavior-preserving: `_write_read_hits_tsv(tmp, ['hit1'])`
  → `_write_extract_tsvs(tmp, ['hit1'], sample_name='')` (already the delegation
  target, writes both TSVs); the empty-targets case still raises via the kept
  guard. `_iter_fastq_read_ids` → `_iter_fastq_records` works since the latter
  raises the same `ValueError`s; adapt the test to consume `(read_id, length)`
  tuples.

### Finding 6 — trivial, fine
- Ready as written.

### Finding 7 — correctly scoped as follow-up
- The cautious framing is right. The barcode test should assert specifically on
  barcode/UMI survival through Picard `SamToFastq` + interleaving
  (`kallisto.py:220-236`): that path emits plain 4-line FASTQ records, so
  barcodes carried in BAM tags (e.g. `CB`/`UB`) will not survive `SamToFastq`
  unless explicitly extracted. Working prior (low confidence): single-cell from
  BAM is in fact broken and `--technology` should be gated to `bulk` for BAM
  input — but confirm with the test before acting.

## Revised readiness

| Finding | Status |
| --- | --- |
| 3, 4, 6 | Ready to implement as written. |
| 5 | Ready, with the public-method caveat noted in the PR. |
| 1 + 2 | Rework required: merge into one input-handling change that (a) hoists a missing-file check covering BAM and FASTQ, (b) DRYs the empty-output headers, and (c) decides + documents the empty-via-derive-targets behavior. As written they ship an incoherent contract. |
| 7 | Remains a tracked follow-up. |
