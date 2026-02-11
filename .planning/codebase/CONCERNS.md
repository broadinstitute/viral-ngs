# Codebase Concerns

**Analysis Date:** 2026-02-11

## Tech Debt

**VCF Record Handling (Overlapping VCF Records):**
- Issue: Overlapping VCF record processing is known to have bugs but lacks fixes
- Files: `src/viral_ngs/phylo/vcf.py` (line 406), `src/viral_ngs/assemble/vcf.py` (line 406)
- Impact: VCF consensus sequence generation may produce incorrect results when VCF records overlap
- Fix approach: Implement proper overlap detection and resolution logic; add comprehensive test coverage for overlapping variants before enabling in production

**VCF Multi-Base Replacement Logic:**
- Issue: Complex case where alleles replace multiple reference bases is untested
- Files: `src/viral_ngs/phylo/vcf.py` (line 415), `src/viral_ngs/assemble/vcf.py` (line 415)
- Impact: Indels and complex variants replacing multiple bases may be incorrectly processed
- Fix approach: Write test cases covering the loop logic at lines 416-426; validate against known problematic VCF inputs

**Deprecated NCBI Tool (tbl2asn):**
- Issue: tbl2asn is deprecated by NCBI but wrapper was removed during monorepo migration; table2asn replacement needs implementation
- Files: `src/viral_ngs/ncbi.py` (line 721-724)
- Impact: GenBank file preparation is skipped with warning; users must manually run table2asn or implement wrapper
- Fix approach: Implement table2asn wrapper matching the different command-line interface; update documentation pointing to manual workaround

**Barcode Handling Gaps:**
- Issue: Missing guardrails for null barcode_3 values and mixed present/absent barcode_3 rows in sample sheets
- Files: `src/viral_ngs/illumina.py` (line 3058)
- Impact: Inline barcode demultiplexing may fail silently or produce incorrect binning when sample sheets have incomplete barcode data
- Fix approach: Add explicit validation and error handling for barcode consistency; document expected sample sheet format

**BWA Alignment I/O Inefficiency:**
- Issue: Read group handling breaks BAM input apart unnecessarily; should use Bwa().mem() directly and pipe to samtools
- Files: `src/viral_ngs/taxon_filter.py` (lines 671-673)
- Impact: Performance degradation and temporary disk space usage; opportunity for optimization
- Fix approach: Refactor to use Bwa().mem() with Popen background process and piped input; validate performance gains

**Barcode Index/Reference Complexity:**
- Issue: Aligner superclass design needed as third aligner added; common functionality duplicated between BWA and Novoalign wrappers
- Files: `src/viral_ngs/core/bwa.py` (line 133-136)
- Impact: Maintenance burden; changes to alignment logic require updates in multiple places
- Fix approach: Extract common aligner interface supporting arbitrary aligners with read group preservation

**H5AD File Handling (Single Sample Limitation):**
- Issue: kb extraction method assumes single-row h5ad (1 sample) and lacks robust multi-sample handling
- Files: `src/viral_ngs/metagenomics.py` (line 1091)
- Impact: Multi-sample h5ad files cannot be processed; scaling to batch processing is blocked
- Fix approach: Redesign extraction to iterate per-sample or process all rows; test with multi-sample h5ad files

**RevertSam Output Handling:**
- Issue: Ambiguous behavior when RevertSam output requested but not needed for unaligned input
- Files: `src/viral_ngs/read_utils.py` (line 256)
- Impact: Silent behavior change (touching empty file) may confuse downstream consumers
- Fix approach: Decide and document: either error out explicitly or validate that empty file is acceptable downstream

**Read Utilities Duplication:**
- Issue: Alignment code duplicated in reports.py instead of using read_utils.py::align_and_fix
- Files: `src/viral_ngs/reports.py` (line 808)
- Impact: Different alignment behavior between reports and main pipeline; maintenance burden
- Fix approach: Consolidate to use align_and_fix; verify reports still produce expected output

**Sample Filtering Logic (Glob Pattern):**
- Issue: Sample name filtering uses fragile dot-based logic depending on filename conventions
- Files: `src/viral_ngs/reports.py` (lines 71-75)
- Impact: Renaming or restructuring sample files can break filtering; edge cases with dots in names not handled
- Fix approach: Replace with explicit whitelist/configuration or proper filename parsing

**KMC Reverse-Complement Handling:**
- Issue: No reverse-complementing of read2 when database is single-strand; unclear if needed
- Files: `src/viral_ngs/classify/kmc.py` (line 249)
- Impact: Potential missed classifications or incorrect k-mer matching for reverse strand
- Fix approach: Document expected behavior; add test with single-strand databases

## Known Bugs

**VCF Parsing Error Handling:**
- Symptoms: VCF parsing errors log location but then re-raise without context
- Files: `src/viral_ngs/intrahost.py` (line 1123-1125), `src/viral_ngs/assembly.py` (line 1277-1279)
- Trigger: Malformed VCF rows with missing CHROM, POS, or FWS fields
- Workaround: Pre-validate VCF files with bcftools before processing

**SnpEff Temp Config File Cleanup:**
- Symptoms: Temporary config file may not be cleaned up if database build fails; exception raised but file remains
- Files: `src/viral_ngs/phylo/snpeff.py` (lines 121-127)
- Trigger: snpEff build command fails after temp config creation
- Workaround: Manually clean /tmp directory; could improve with context manager

## Security Considerations

**File Handle Leaks:**
- Risk: Multiple places use open() without context managers; files may not be closed on exceptions
- Files:
  - `src/viral_ngs/file_utils.py`: csv.DictReader with open() without context manager
  - `src/viral_ngs/core/samtools.py`: Lines ~456, ~457 (open for stdout/stderr without context)
  - `src/viral_ngs/core/sambamba.py`: open() without context manager for file handles
  - `src/viral_ngs/core/minimap2.py`: open() without context manager
- Current mitigation: Process cleanup on normal completion; file handles eventually close
- Recommendations: Wrap all open() calls in with statements or use util_file helpers; audit for resource leaks

**Bare Except Clauses:**
- Risk: Bare except: clauses catch and suppress KeyboardInterrupt and SystemExit
- Files:
  - `src/viral_ngs/intrahost.py` (line 1123)
  - `src/viral_ngs/phylo/snpeff.py` (line 124)
  - `src/viral_ngs/assembly.py` (line 1277)
- Current mitigation: Most bare excepts are followed by re-raise or logging
- Recommendations: Replace all bare excepts with Exception or specific exception types; add linting rule to prevent

**Path Construction Safety:**
- Risk: Some path concatenation uses string operations rather than os.path.join consistently
- Files: `src/viral_ngs/reports.py` (lines 71-76) uses string concatenation with sample + ".*.bam"
- Current mitigation: Paths validated by glob before use
- Recommendations: Audit all path construction; prefer pathlib.Path for new code

## Performance Bottlenecks

**Illumina Module Size (3826 lines):**
- Problem: Single file handles demultiplexing, sample sheet parsing, barcode assignment, and reporting
- Files: `src/viral_ngs/illumina.py`
- Cause: Growth over time without refactoring into separate modules
- Improvement path: Split into: illumina_demux.py, illumina_samplesheet.py, illumina_barcodes.py; establish clear interfaces

**Illumina Indices Lookup Table (2742 lines):**
- Problem: Large hardcoded index reference data; lookup performance not profiled
- Files: `src/viral_ngs/core/illumina_indices.py`
- Cause: Comprehensive enumeration of Illumina indices as inline dictionaries
- Improvement path: Consider loading from external data file or generating at build time; profile dictionary lookups

**Core Misc Utilities (1361 lines):**
- Problem: Miscellaneous utilities accumulating without clear module boundaries
- Files: `src/viral_ngs/core/misc.py`
- Cause: Catch-all module for functions without specific homes
- Improvement path: Categorize functions (file handling, subprocess, string manipulation); move to specific modules

**Feature Table Types (1297 lines):**
- Problem: Large feature table type definitions hard to modify; syntax complex
- Files: `src/viral_ngs/phylo/feature_table_types.py`
- Cause: Comprehensive GenBank feature type reference
- Improvement path: Consider data-driven approach with external definition file

**File Utility Module (1257 lines):**
- Problem: Mixed concerns: path operations, compression, archive handling, I/O utilities
- Files: `src/viral_ngs/core/file.py`
- Cause: Centralization to avoid duplication
- Improvement path: Extract compression handling to separate module; ensure single responsibility

## Fragile Areas

**Illumina Demultiplexing Pipeline:**
- Files: `src/viral_ngs/illumina.py` (full module), `src/viral_ngs/core/illumina_utils.py`
- Why fragile: Complex state machines for barcode assignment; multiple interdependent data structures; XML parsing of RunInfo
- Safe modification: Use existing test fixtures in tests/unit/core/test_illumina.py; validate with real RunInfo examples
- Test coverage: demux logic tested but barcode assignment edge cases need expansion

**VCF Processing:**
- Files: `src/viral_ngs/phylo/vcf.py`, `src/viral_ngs/assemble/vcf.py` (duplicated)
- Why fragile: Complex string manipulation; overlapping record handling incomplete; untested multi-base cases
- Safe modification: Add unit tests before changes; validate output against bcftools consensus
- Test coverage: Basic functionality tested; edge cases (overlaps, indels, low-complexity) not covered

**Subprocess Management:**
- Files: `src/viral_ngs/core/misc.py` (subprocess utilities)
- Why fragile: Multiple Popen calls with complex piping; stdout/stderr handling inconsistent
- Safe modification: Use existing wrapper functions; test error cases with failing subprocess
- Test coverage: Happy path tested; subprocess failure scenarios incomplete

**Read Group Preservation in Aligners:**
- Files: `src/viral_ngs/core/bwa.py`, `src/viral_ngs/core/novoalign.py`
- Why fragile: Picard tools coordinate with samtools for read group handling; complex file conversions
- Safe modification: Add comprehensive test with multi-group BAM input
- Test coverage: Single read group tested; multi-group and missing @RG scenarios untested

**SNP Effect Prediction:**
- Files: `src/viral_ngs/phylo/snpeff.py`
- Why fragile: Temporary directory and config file management; database download dependencies
- Safe modification: Add integration test with known SNPs; validate output format
- Test coverage: Basic execution tested; database edge cases and missing databases not covered

## Scaling Limits

**Multi-Sample H5AD Processing:**
- Current capacity: Single-sample h5ad files only (due to assumption at line 1091 of metagenomics.py)
- Limit: Cannot process batches of samples; each sample requires separate run
- Scaling path: Refactor extraction to handle multiple rows; batch processing at higher level

**Concurrent File Handle Limits:**
- Current capacity: Using ProcessPoolExecutor with max_workers=sanitize_thread_count()
- Limit: Not profiled; potential exhaustion with very large concurrent jobs
- Scaling path: Profile with large thread counts; consider resource pooling or queueing

**Memory Usage in Large Arrays:**
- Current capacity: Illumina indices and feature tables loaded entirely into memory
- Limit: Unknown; not profiled with large reference datasets
- Scaling path: Monitor memory usage; consider lazy-loading or external data stores

**Temporary Disk Space:**
- Current capacity: Multiple intermediate BAM files created during demultiplexing and alignment
- Limit: Pipeline may fail if /tmp fills up; no cleanup on failure
- Scaling path: Profile disk usage; implement streaming where possible; ensure cleanup on error

## Dependencies at Risk

**tbl2asn Deprecation:**
- Risk: NCBI deprecated tbl2asn in favor of table2asn; old wrapper removed
- Impact: GenBank submission workflow broken
- Migration plan: Implement table2asn wrapper with new interface; document manual workaround

**Pysam Version Compatibility:**
- Risk: Old comments suggest pysam 0.8.1 had bugs with nosetests; code predates pysam 0.9.1 improvements
- Impact: Potential performance issues if using old pysam; testing framework changed (pytest, not nosetests)
- Migration plan: Verify minimum pysam version in requirements; consider consolidating samtools usage

**Arrow Parser for Dates:**
- Risk: arrow.parser.ParserError used in two places; library maintenance status should be verified
- Impact: Date parsing failures in RunInfo or sample sheet parsing
- Migration plan: Consider Python datetime library if possible; add error handling

**Matplotlib in Illumina Module:**
- Risk: Matplotlib imported in non-GUI context; may fail in headless environments
- Impact: Import errors if matplotlib backend not configured
- Migration plan: Make matplotlib import lazy or conditional; consider alternative plotting libraries

## Missing Critical Features

**H5AD Batch Processing:**
- Problem: Cannot process multiple samples from kb extraction in single command
- Blocks: Scaling metagenomics classification to large sample batches
- Workaround: Run separately per sample; parallel via external orchestration

**IUPAC Ambiguity Code Handling:**
- Problem: Heterozygous calls not converted to IUPAC codes; marked TODO but not implemented
- Files: `src/viral_ngs/phylo/vcf.py` (line 400-401)
- Blocks: Consensus sequences may show arbitrary genotype instead of ambiguity
- Workaround: Post-process VCF or consensus with external tool

**Better Sample Name Filtering:**
- Problem: Glob-based filtering fragile and undocumented
- Files: `src/viral_ngs/reports.py` (lines 71-75)
- Blocks: Reports generation can fail with unusual sample naming
- Workaround: Rename samples to follow expected pattern (no dots)

**R2 Barcode Support:**
- Problem: R2 barcode support commented out in splitcode; not implemented
- Files: `src/viral_ngs/core/splitcode.py` (line 909)
- Blocks: Cannot demux using reverse-read barcodes
- Workaround: Use external demultiplexing tool

## Test Coverage Gaps

**VCF Consensus Generation:**
- What's not tested: Overlapping VCF records, multi-base indel handling, complex variant scenarios
- Files: `src/viral_ngs/phylo/vcf.py`, `src/viral_ngs/assemble/vcf.py`
- Risk: Untested code paths may silently produce incorrect consensus sequences
- Priority: High

**Illumina Barcode Assignment Edge Cases:**
- What's not tested: Missing barcode_3, mixed present/absent barcode_3 rows, barcode collisions
- Files: `src/viral_ngs/illumina.py` (lines 3000-3070)
- Risk: Edge cases may cause crashes or incorrect sample binning
- Priority: High

**Subprocess Failure Handling:**
- What's not tested: Tool execution failures, missing executables, permission errors
- Files: `src/viral_ngs/core/misc.py` (subprocess wrappers)
- Risk: Errors not propagated clearly; potential for silent failures
- Priority: Medium

**File Resource Cleanup:**
- What's not tested: Exceptions during file operations; cleanup of temporary files
- Files: Multiple (file_utils.py, samtools.py, sambamba.py, minimap2.py)
- Risk: Disk space leaks; file handle leaks under error conditions
- Priority: Medium

**Read Group Preservation Multi-Group:**
- What's not tested: Multiple read groups in single BAM; complex RG tag scenarios
- Files: `src/viral_ngs/core/bwa.py`, `src/viral_ngs/core/novoalign.py`
- Risk: Alignment may fail or lose read group information with multi-group input
- Priority: Medium

**SNPEff Database Edge Cases:**
- What's not tested: Missing databases, corrupted database files, download failures
- Files: `src/viral_ngs/phylo/snpeff.py`
- Risk: Silent annotation failures or confusing error messages
- Priority: Low

**Multi-Aligner Refactoring:**
- What's not tested: Third aligner integration; interface stability
- Files: `src/viral_ngs/core/bwa.py`, `src/viral_ngs/core/novoalign.py`
- Risk: Duplicate code maintenance burden; inconsistent behavior across aligners
- Priority: Low

---

*Concerns audit: 2026-02-11*
