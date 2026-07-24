Lyra classifier post-processing
===============================

This page is the normative pipeline-author contract for converting native
Lyra scores into normalized fragment classifications, a run summary, and a
viral-read BAM.

Runtime boundary
----------------

viral-ngs consumes Lyra output only. Lyra inference remains in the separate
``lyra-cuda`` image, where external workflow code invokes ``lyra_cli``.
viral-ngs neither executes nor packages Lyra, PyTorch, CUDA, model code or
checkpoints, or a GPU runtime. It does not contain the viral-pipelines WDL
implementation that invokes inference.

Command
-------

The default-threshold invocation is:

.. code-block:: text

   metagenomics lyra_postprocess read_scores.tsv source.bam SAMPLE_ID normalized.tsv summary.tsv viral.bam

Omitting ``--threshold`` uses the inclusive threshold ``0.8``. A compressed
invocation with an explicit threshold and scratch base is:

.. code-block:: text

   metagenomics lyra_postprocess read_scores.tsv.zst source.bam SAMPLE_ID normalized.tsv.gz summary.tsv viral.bam --threshold 0.95 --tmp_dir /path/to/scratch

All six positional values are required in this order: native score table,
source BAM, sample identifier, normalized TSV, summary TSV, and viral BAM.
``--threshold`` is the only classification-threshold option; it must be a
finite value from 0 through 1, inclusive. Native score tables and normalized
tables support exact lowercase ``.tsv``, ``.tsv.gz``, and ``.tsv.zst``
suffixes. The summary path must end in exact lowercase ``.tsv``, and the viral
output path must end in exact lowercase ``.bam``. A successful command writes
nothing to stdout. ``--tmp_dir`` selects the base for command-managed SQLite
scratch state.

Native score input
------------------

The native table must have exactly these ordered columns:

.. list-table:: Native Lyra score schema
   :header-rows: 1
   :widths: 10 25 65

   * - Order
     - Column
     - Meaning
   * - 1
     - ``read_id``
     - The exact, unmodified BAM QNAME for one Lyra-scored BAM record.
   * - 2
     - ``score``
     - A strict ordinary finite decimal from 0 through 1, inclusive.
   * - 3
     - ``call``
     - The literal native value ``0`` or ``1``; it is validated but never
       trusted when computing output calls.

One data row corresponds to one BAM record scored by Lyra. Malformed text,
headers, row widths, identifiers, scores, and native calls are rejected rather
than repaired or silently accepted. QNAME identity is exact and
case-sensitive: identifiers are not trimmed, suffix-stripped, whitespace-split,
case-normalized, or otherwise modified.

BAM reconciliation
------------------

A BAM record is Lyra-eligible exactly when it is primary, non-supplementary,
has a stored query sequence, and that sequence is at least 50 bases long.
For every exact QNAME, score-row multiplicity must equal eligible-BAM-record
multiplicity. Eligible paired evidence must also have coherent, exclusive mate
roles: a complete pair has exactly one read-1-only record and one read-2-only
record, while one eligible paired record with exactly one mate role is an
incomplete pair. Ineligible records do not contribute score or pairing
evidence.

Fragment calling
----------------

Calls are recomputed only from the validated scores and configured inclusive
threshold:

.. list-table:: Fragment calling truth table
   :header-rows: 1
   :widths: 25 45 30

   * - Pairing
     - Score condition
     - Call
   * - ``Single-end``
     - Its one score is ``>= threshold``
     - ``Viral``
   * - ``Paired-complete``
     - Both scores are ``>= threshold``
     - ``Viral``
   * - ``Paired-incomplete``
     - Any available score
     - ``Non-viral``

All other score conditions for single-end and complete-pair fragments produce
``Non-viral``. Equality with the threshold passes. Fragment calls are
independent of score-row and BAM-record order. Final normalized rows are sorted
by exact UTF-8 ordinal ``READ_ID``.

Normalized classifications
--------------------------

The normalized table has exactly these ordered columns:

.. list-table:: Normalized classification schema
   :header-rows: 1
   :widths: 10 30 60

   * - Order
     - Column
     - Meaning
   * - 1
     - ``SAMPLE_ID``
     - The validated sample identifier supplied to the command.
   * - 2
     - ``READ_ID``
     - The exact, unmodified fragment QNAME.
   * - 3
     - ``LYRA_N_SCORES``
     - Number of eligible Lyra score records for the fragment: one or two.
   * - 4
     - ``LYRA_PAIRING``
     - Reconciled pairing state.
   * - 5
     - ``LYRA_MIN_SCORE``
     - Minimum validated score for the fragment.
   * - 6
     - ``LYRA_MAX_SCORE``
     - Maximum validated score for the fragment.
   * - 7
     - ``LYRA_THRESHOLD``
     - Inclusive threshold applied when recomputing the call.
   * - 8
     - ``LYRA_CALL``
     - Recomputed fragment-level call.

The ``LYRA_PAIRING`` vocabulary is exactly ``Single-end``,
``Paired-complete``, and ``Paired-incomplete``. The ``LYRA_CALL`` vocabulary is
exactly ``Viral`` and ``Non-viral``. Exactly one row is emitted for every
``(SAMPLE_ID, READ_ID)``, including every nonviral fragment. Scores and the
threshold use canonical ordinary decimal text without exponent notation or
redundant trailing fractional zeroes. Plain, gzip, and Zstandard normalized
forms have identical logical UTF-8 content with LF record terminators.

Summary
-------

The uncompressed summary contains exactly one data row under these ordered
columns:

.. list-table:: Summary schema
   :header-rows: 1
   :widths: 10 35 55

   * - Order
     - Column
     - Meaning
   * - 1
     - ``SAMPLE_ID``
     - The validated sample identifier supplied to the command.
   * - 2
     - ``LYRA_THRESHOLD``
     - The inclusive threshold applied to fragment scores.
   * - 3
     - ``LYRA_INPUT_BAM_RECORDS``
     - Number of all records read from the source BAM.
   * - 4
     - ``LYRA_ELIGIBLE_BAM_RECORDS``
     - Number of source records eligible for Lyra score reconciliation.
   * - 5
     - ``LYRA_SCORE_RECORDS``
     - Number of validated native score rows.
   * - 6
     - ``LYRA_FRAGMENTS``
     - Number of reconciled exact-QNAME fragments and normalized data rows.
   * - 7
     - ``LYRA_SINGLE_END_FRAGMENTS``
     - Number of fragments classified as ``Single-end``.
   * - 8
     - ``LYRA_COMPLETE_PAIR_FRAGMENTS``
     - Number of fragments classified as ``Paired-complete``.
   * - 9
     - ``LYRA_INCOMPLETE_PAIR_FRAGMENTS``
     - Number of fragments classified as ``Paired-incomplete``.
   * - 10
     - ``LYRA_VIRAL_FRAGMENT_CALLS``
     - Number of fragments called ``Viral``.
   * - 11
     - ``LYRA_NONVIRAL_FRAGMENT_CALLS``
     - Number of fragments called ``Non-viral``.
   * - 12
     - ``LYRA_OUTPUT_BAM_RECORDS``
     - Actual number of records emitted to the viral BAM.

Every successful summary satisfies all of these cardinality equations:

.. code-block:: text

   LYRA_ELIGIBLE_BAM_RECORDS = LYRA_SCORE_RECORDS
   LYRA_SCORE_RECORDS = LYRA_SINGLE_END_FRAGMENTS + LYRA_INCOMPLETE_PAIR_FRAGMENTS + 2 * LYRA_COMPLETE_PAIR_FRAGMENTS
   LYRA_FRAGMENTS = LYRA_SINGLE_END_FRAGMENTS + LYRA_COMPLETE_PAIR_FRAGMENTS + LYRA_INCOMPLETE_PAIR_FRAGMENTS
   LYRA_FRAGMENTS = LYRA_VIRAL_FRAGMENT_CALLS + LYRA_NONVIRAL_FRAGMENT_CALLS
   0 <= LYRA_OUTPUT_BAM_RECORDS <= LYRA_INPUT_BAM_RECORDS
   LYRA_OUTPUT_BAM_RECORDS = 0 exactly when LYRA_VIRAL_FRAGMENT_CALLS = 0

Viral BAM
---------

The viral BAM contains every source record whose exact QNAME belongs to a
``Viral`` fragment. This includes same-QNAME ineligible, secondary,
supplementary, short, or sequence-less companion records that did not
contribute calling evidence. Filtering preserves the source header and record
order, and preserves flags, tags, sequence, qualities, and other record fields.
The command performs no sorting, indexing, tagging, or addition of a program
record.

Empty and no-hit runs
---------------------

A valid empty run emits a normalized header with no data rows, one complete
summary row with zero counts, and a header-preserving zero-record BAM. For an
all-ineligible nonempty source BAM, ``LYRA_INPUT_BAM_RECORDS`` can be nonzero
while the eligible, score, fragment, call, and output counts remain zero. A
valid no-hit run retains every normalized fragment row as ``Non-viral`` and
emits a header-preserving zero-record BAM.

Publication and failure behavior
--------------------------------

Pre-existing final entries, input/output aliases, and output/output aliases are
rejected without clobbering caller-owned paths. The command writes hidden
same-parent stages and performs format-aware full readback validation and, on
supported Linux filesystems, durability barriers before publication.
Normalized output publishes first, the viral BAM second, and the summary last
as the completion marker.

On any reported failure, the process's nonzero status is authoritative. For a
catchable failure, successful cleanup syscalls remove run-owned final entries
and hidden stages. If rollback unlink fails, the observed residual is preserved
and bounded, ordered diagnostics report it without replacing the authoritative
primary cause. ``SIGKILL``, host loss, and power loss can also leave residual
entries. Therefore summary-path presence after a reported failed rollback is
not success; the summary is a successful completion marker only together with
a zero process status.
