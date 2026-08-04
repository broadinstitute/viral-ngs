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
- When new fixable HIGH/CRITICAL CVEs are detected (i.e. CVE IDs not already present in any open or closed GH issue title), the workflow invokes Claude Sonnet 5 on Vertex AI to triage each one and files a GitHub issue per CVE (labels: `security`, `cve`)
- **Auto-fix PR (`cve-fix-pr` job)**: after new issues are filed, Claude assesses *all* currently-open `cve` issues and either opens a single `cve-fix` + `security` PR — only if every open issue meets the moderate confidence bar (version-floor bump, precedent-mirroring Dockerfile mitigation, or precedent-mirroring `.trivyignore` entry; all-or-nothing) — or uploads a `SKIPPED.md` decision-log artifact. The PR is authored with a GitHub App token so it triggers the required CI checks.
- **Recovery / manual re-run**: the auto-fix job only fires when the scan finds a *new* CVE, so once an issue is filed the CVE is no longer "new" and a failed or skipped fix run will not retry on its own. To re-run it on demand against all currently-open `cve` issues, dispatch manually: `gh workflow run container-scan.yml -f force_fix_pr=true`
- **Conda availability probe**: both Claude steps are preceded by `.github/scripts/probe-conda-availability.py`, which queries the conda-forge and bioconda APIs and writes `.cve-fix-context/conda-availability.json`. Per affected package it reports whether the Trivy `FixedVersion` exists on a channel, whether it is built for **both** `linux-64` and `linux-aarch64` (or is `noarch`), whether a requirements pin can fix the finding at all (`fix_shape` — an Ubuntu `os-package` or a `bundled-jar` cannot), and which other pinned packages carry a constraint that would exclude the fix version (`constrained_by`, e.g. `pyopenssl<50` capping `cryptography`). Uploaded as a workflow artifact. Before this existed the agents had no network tool at all and skipped *every* version-floor bump as unverifiable — the primary category the pipeline was built for.
- **CI verifies ARM64, so the agent does not have to**: `docker.yml` builds all six tiers for both amd64 and native arm64 on every PR and re-scans each with Trivy, so a solver conflict from a floor bump fails required checks before merge. The auto-fix prompt says so explicitly and instructs the agent to reserve SKIP for "I cannot determine *what* the fix is", not "I cannot prove the fix works". Unproven solvability used to be the standard reason for skipping.
- **Steering the bot via issue comments**: maintainer comments on a `cve` issue are passed to the auto-fix agent and outrank the machine-generated issue body — use them to specify or amend a fix. Comments are hard-filtered to `OWNER`/`MEMBER`/`COLLABORATOR` in the workflow's `jq`, not in the prompt: this is a public repo, so an unfiltered comments array would be an untrusted-instruction channel into an agent that edits files and opens PRs.
- **Agent guardrails**: the auto-fix agent has domain-scoped `WebFetch` (no `WebSearch`, since its output is committed); the staged diff is checked against a path allowlist (`docker/requirements/*.txt`, `docker/Dockerfile.*`, `.trivyignore`) before commit; and the GitHub App token is kept out of `.git/config` while the agent runs (`persist-credentials: false`, push remote attached afterwards). The triage agent gets `WebSearch` too — it cannot edit or push.
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
3. **Check if a fix version exists on conda-forge/bioconda** -- not just on PyPI. Run
   `.github/scripts/probe-conda-availability.py --trivy-json trivy-results.json --cve-ids "<CVE-ID>" --out /tmp/avail.json`
   (works locally against a downloaded scan artifact) rather than checking by hand; it also
   reports which other pinned packages constrain the version you want.
4. **Test on ARM64** -- solver conflicts are the most common failure mode. Opening a PR is
   a legitimate way to test this: `docker.yml` builds native arm64 for every tier.
5. **If the fix version conflicts**: consider whether the CVE is exploitable in your deployment model. Document the risk assessment in `.trivyignore` or `vulnerability-mitigation-status.md`.
6. **If the vulnerable code is unused**: delete the binary/file inline in the Dockerfile (same RUN layer as install to avoid bloating images)

## Key Files

| File | Purpose |
|------|---------|
| `.trivy-ignore-policy.rego` | Rego policy for class-level CVE filtering |
| `.trivyignore` | Per-CVE exceptions with justifications |
| `.github/workflows/docker.yml` | Build-time scanning (SARIF + JSON), amd64 + native arm64 |
| `.github/workflows/container-scan.yml` | Daily scheduled scanning, Claude triage, auto-fix PR |
| `.github/scripts/probe-conda-availability.py` | Pre-agent conda-forge/bioconda availability + reverse-constraint probe |
| `vulnerability-mitigation-status.md` | Local-only tracking doc (not committed) |
