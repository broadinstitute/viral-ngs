# Container Vulnerability Management

Guidance for scanning, triaging, and mitigating container image vulnerabilities
in the viral-ngs Docker image hierarchy.

## Scanning

Container images are scanned for vulnerabilities using [Trivy](https://aquasecurity.github.io/trivy/):

### PR/Push Scans (`docker.yml`)

- Runs on every PR/push after image build
- **SARIF scan**: CRITICAL/HIGH, fixable-only → uploads to GitHub Security tab
- **JSON scan**: CRITICAL/HIGH, fixable-only, **exit-code: 1** → gates CI (fails if CVEs found)
- Both apply `.trivyignore` and `.trivy-ignore-policy.rego`

### Scheduled Scans (`container-scan.yml`)

- Runs daily at 09:00 UTC on published `main-mega-amd64` image
- **SARIF scan**: CRITICAL/HIGH, fixable-only → uploads to GitHub Security tab
- **JSON scan**: CRITICAL/HIGH/MEDIUM, all vulnerabilities (fixable + non-fixable), **exit-code: 0** → uploaded as artifact for comprehensive debugging and Claude analysis
- When new fixable HIGH/CRITICAL CVEs are detected (i.e. CVE IDs not already present in any open or closed GH issue title), the workflow invokes Claude Sonnet 4.6 on Vertex AI to triage each one and files a GitHub issue per CVE (labels: `security`, `cve`)
- The Vertex/WIF infra used here is documented in `.agents/skills/claude-on-vertex-ci/SKILL.md`

**Key difference**: The scheduled scan's JSON output includes MEDIUM severity and non-fixable vulnerabilities to provide full visibility for risk assessment and debugging, while the control flow (issue filing, CI gates) still filters to fixable HIGH/CRITICAL only.

Per-CVE exceptions go in `.trivyignore` with mandatory justification comments.

## Rego Policy (`.trivy-ignore-policy.rego`)

The Rego policy filters CVEs that are architecturally inapplicable to ephemeral batch containers:

- **AV:P** (Physical access required) -- containers are cloud-hosted
- **AV:A** (Adjacent network required) -- no attacker on same network segment
- **AV:L + UI:R** (Local + user interaction) -- no interactive sessions
- **AV:L + PR:H** (Local + high privileges) -- containers run non-root
- **AV:L + S:U** (Local + scope unchanged) -- attacker already has code execution and impact stays within the ephemeral container

Changes to this policy should be reviewed carefully. The comments in the file explain the rationale and risk for each rule.

## Common Vulnerability Sources

**Python transitive deps**: Pin minimum versions in `docker/requirements/*.txt`. Prefer conda packages over pip. Check conda-forge availability before assuming a version exists -- conda-forge often lags PyPI by days/weeks.

**Java fat JARs** (picard, gatk, snpeff, fgbio): Bioinformatics Java tools are distributed as uber JARs with all dependencies bundled inside. Trivy detects vulnerable libraries (log4j, commons-compress, etc.) baked into these JARs. Version bumps can cause ARM64 conda solver conflicts because Java tools pull in openjdk -> harfbuzz -> icu version chains that clash with other packages (r-base, boost-cpp, pyicu). Always check:
1. Whether the tool is actually flagged by Trivy (don't bump versions unnecessarily)
2. Whether the CVE applies (e.g., log4j 1.x is NOT vulnerable to Log4Shell)
3. Whether the desired version resolves on ARM64 before pushing

**Go binaries**: Some conda packages bundle compiled Go binaries (e.g., mafft's `dash_client`, google-cloud-sdk's `gcloud-crc32c`). If the binary is unused, delete it in the Dockerfile. Delete from **both** the installed location and `/opt/conda/pkgs/*/` (conda package cache) -- Trivy scans the full filesystem.

**Vendored copies**: Packages like google-cloud-sdk and setuptools bundle their own copies of Python libraries that may be older than what's in the conda environment. Trivy flags these vendored copies separately. Options: delete the vendored directory (if not needed at runtime), or accept the risk in `.trivyignore` with justification.

## ARM64 Solver Conflicts

The conda solver on ARM64 (linux-aarch64) is more constrained than amd64 because fewer package builds exist. Common conflict patterns:

- **icu version conflicts**: Many packages (openjdk, r-base, boost-cpp, pyicu) pin specific icu version ranges. Bumping one package can make the entire environment unsolvable.
- **libdeflate/htslib conflicts**: lofreq 2.1.5 pins old htslib/libdeflate versions that conflict with newer pillow/libtiff.
- **openjdk version escalation**: snpeff 5.2+ requires openjdk>=11, 5.3+ requires openjdk>=21. Higher openjdk versions pull in harfbuzz->icu chains that conflict with everything.

When a solver conflict occurs: revert the change, check what version the solver was picking before, and pin to that exact version if it already addresses the CVE.

## Mitigation Decision Process

When triaging a CVE:

1. **Check the CVSS vector** -- does the Rego policy already filter it?
2. **Identify the source package** -- use Trivy JSON output (`PkgName`, `PkgPath`, `InstalledVersion`)
3. **Check if a fix version exists on conda-forge/bioconda** -- not just on PyPI
4. **Test on ARM64** -- solver conflicts are the most common failure mode
5. **If the fix version conflicts**: consider whether the CVE is exploitable in your deployment model. Document the risk assessment in `.trivyignore` or `vulnerability-mitigation-status.md`.
6. **If the vulnerable code is unused**: delete the binary/file inline in the Dockerfile (same RUN layer as install to avoid bloating images)

## Key Files

| File | Purpose |
|------|---------|
| `.trivy-ignore-policy.rego` | Rego policy for class-level CVE filtering |
| `.trivyignore` | Per-CVE exceptions with justifications |
| `.github/workflows/docker.yml` | Build-time scanning (SARIF + JSON) |
| `.github/workflows/container-scan.yml` | Weekly scheduled scanning |
| `vulnerability-mitigation-status.md` | Local-only tracking doc (not committed) |
