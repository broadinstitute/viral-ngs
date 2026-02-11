# External Integrations

**Analysis Date:** 2026-02-11

## APIs & External Services

**NCBI/GenBank:**
- Service: NCBI GenBank and Genbank CoreNucleotide database
- What it's used for: Sequence retrieval, taxonomy lookups, accession validation
  - SDK/Client: BioPython `Bio.Entrez` module
  - Location: `src/viral_ngs/phylo/genbank.py` (function `_fetch_from_nuccore()`)
  - Auth: Email address (required by NCBI policy for rate limiting)
  - Optional: `api_key` parameter for faster requests
  - Usage: `Entrez.email = emailAddress`, `Entrez.api_key = api_key` if provided

**NCBI eSRA (Sequence Read Archive):**
- Service: NCBI SRA for read archival
- What it's used for: Submitting sequencing data to public archives
  - Location: `src/viral_ngs/ncbi.py` (module contains NCBI submission utilities)
  - Auth: Credentials managed externally (not embedded in code)

**NCBI Taxonomy Database:**
- Service: NCBI Taxonomy dump files (downloadable)
- What it's used for: Taxonomic classification and filtering
  - Location: `src/viral_ngs/classify/taxonomy.py`
  - Input: Pre-downloaded files:
    - `gi->taxid` mapping (gi.dmp)
    - `names.dmp` (taxonomy names)
    - `nodes.dmp` (taxonomy hierarchy)
    - `taxdump.tar.gz` (complete taxonomy dump)
  - Usage: Command-line arguments `--taxdump_tar_gz` and optional `--get_accessions` flag
  - Note: Taxonomy files are loaded locally, not fetched at runtime

## Data Storage

**Databases:**
- Type: No persistent database system
- File format storage: Primary use of file-based data
  - Sequence formats: FASTA, FASTQ, SAM/BAM, VCF/BCF
  - Metadata: YAML, CSV, JSON
  - Compressed: gzip (bgzip), lz4, zstandard (zstd)
- SQLite: Minimal usage - `sqlite3` module imported in `src/viral_ngs/core/file.py` but used only for temporary operations
- Connection: None - file-based only, no network database connections

**File Storage:**
- Local filesystem only - No cloud object storage integration in core code
- Path handling: Standard Python `os.path` module
- Temporary files: `tempfile` module with optional debugging via `VIRAL_NGS_TMP_DIRKEEP` environment variable
- Mounted volumes in Docker for input/output data

**Caching:**
- No explicit caching layer (in-memory only during execution)
- Tool outputs cached as files on filesystem
- No Redis, Memcached, or similar services

## Authentication & Identity

**Auth Provider:**
- Type: No centralized authentication system
- Approach: Tool-based and environment-specific
  - NCBI Entrez: Email + optional API key (command-line arguments)
  - Novoalign: License file path via environment variable `NOVOALIGN_LICENSE_PATH`
  - AWS/Cloud: Delegated to host credentials (via AWS CLI/SDK in baseimage)
- No user authentication or authorization within the application

## Cloud & Infrastructure

**AWS Integration:**
- Service: AWS Command-line tools available in containers
- Package: awscli 1.32.0+ (installed in baseimage.txt)
- Usage: Optional - for users needing S3 access to input/output files
- Integration: Delegated to AWS credentials on host or container
- No embedded AWS SDK calls in application code

**Google Cloud Integration:**
- Service: Google Cloud Storage client
- Package: google-cloud-storage 2.14.0+ (installed in baseimage.txt)
- Usage: Optional - for users with GCS buckets
- Integration: Delegated to Google Cloud credentials on host
- No embedded GCS calls in application code

**WDL Workflow Execution:**
- Service: miniwdl 1.11.1+ (installed in baseimage)
- Backend: udocker (unprivileged container runtime)
- Configuration: `MINIWDL__SCHEDULER__CONTAINER_BACKEND=udocker` in baseimage
- Usage: For orchestrating viral-ngs tasks in WDL workflows
- Related: Higher-level pipelines at https://github.com/broadinstitute/viral-pipelines

## Monitoring & Observability

**Error Tracking:**
- Type: None (no Sentry, DataDog, or similar integration)
- Approach: Standard Python logging

**Logs:**
- Framework: Python `logging` module (standard library)
- Configuration: CLI-based log level setup via `setup_logger()` in `src/viral_ngs/core/cmd.py`
- Format: `%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s`
- Output: StreamHandler (stderr/stdout)
- Levels: Configurable via command-line `--loglevel` argument (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- No external logging service or structured logging framework

**Metrics/Coverage:**
- codecov.io integration via GitHub Actions
- Configuration: `.codecov.yml` for CI/CD upload settings
- Coverage tracking: `.coveragerc` for pytest-cov configuration

## CI/CD & Deployment

**Container Registry:**
- Primary: Quay.io (`quay.io/broadinstitute/viral-ngs`)
- Secondary: GitHub Container Registry (`ghcr.io/broadinstitute/viral-ngs`)
- Image naming convention: `{main|version}-{baseimage|core|assemble|classify|phylo|mega}`

**CI Pipeline:**
- Service: GitHub Actions
- Workflow file: `.github/workflows/docker.yml`
- Triggers: Push to main branch, all feature branches, version tags, pull requests to main
- Jobs:
  - Path filtering: Smart detection of changed modules (dorny/paths-filter@v3)
  - Version detection: Git describe for semantic versioning
  - Multi-architecture builds: linux/amd64 and linux/arm64
  - Test execution: Per-module pytest runs within Docker containers
  - Registry push: Conditional push to Quay.io and GHCR based on branch/tag
- Version scheme: From git tags (with or without 'v' prefix) or git describe fallback

**Documentation:**
- Service: Read the Docs
- Configuration: `.readthedocs.yml`
- Build tool: Sphinx
- Source: `docs/` directory
- Python build: Ubuntu 22.04, Python 3.10
- Auto-deployment: On push to main branch

## Environment Configuration

**Required environment variables:**
- None for basic operation (all tools accessed via command-line)
- Optional for features:
  - `VIRAL_NGS_SOURCE_DIR` - Override source directory location
  - `NOVOALIGN_LICENSE_PATH` - Path to novoalign license
  - `NCBI_EMAIL` - Email for NCBI queries (passed as argument instead)
  - `NCBI_API_KEY` - API key for NCBI (passed as argument instead)

**Secrets location:**
- Not embedded in code
- Managed externally:
  - AWS credentials: IAM roles or credential files (host-managed)
  - Google Cloud credentials: Service account keys (host-managed)
  - NCBI API keys: Command-line arguments or environment
  - Novoalign licenses: File paths (user-provided)

## Webhooks & Callbacks

**Incoming:**
- None - No webhook endpoints defined
- Tool outputs: Generated as files, not streamed to external services

**Outgoing:**
- None - No outbound callbacks or notifications to external services
- Results are written to local/mounted filesystems

## Network Configuration

**HTTP Calls:**
- NCBI Entrez API: `Bio.Entrez` handles HTTP requests to NCBI servers
  - Endpoint: NCBI E-utilities (e.g., eefetch, esearch)
  - Rate limiting: Email-based throttling (configurable in Entrez object)
  - Chunking: Accession lists chunked to max 500 per request (NCBI guideline compliance)

**DNS/Resolution:**
- Standard Python `urllib.request` for potential file downloads
- No custom DNS configuration needed

---

*Integration audit: 2026-02-11*
