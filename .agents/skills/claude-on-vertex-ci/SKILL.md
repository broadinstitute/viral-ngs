# Claude on Vertex AI from GitHub Actions

Reusable infrastructure for invoking Claude (via [`anthropics/claude-code-action`](https://github.com/anthropics/claude-code-action))
on Google Vertex AI from GitHub Actions workflows in this repo. Authentication
is via Workload Identity Federation — no long-lived secrets.

## When to Use

Add a new Claude-in-CI use case (PR review, CVE triage, dependency-update
analysis, automated docs generation, etc.) when:

- The task is well-suited to an LLM with tool use (reading files, running CLI
  commands, calling `gh`)
- The task fires on a GitHub event (schedule, workflow_dispatch, push, PR)
- You want an issue, PR comment, or artifact as the output

The first use case (and template) is the **CVE triage step in
`container-scan.yml`** — see `.agents/skills/container-vulns/SKILL.md`.

## What's Already Provisioned

**GCP project: `viral-seq-ai`**

| Resource | Identifier |
|---|---|
| Workload Identity Pool | `github-actions-pool` (global) |
| OIDC provider | `broadinstitute-github` |
| Provider attribute condition | `assertion.repository_owner == 'broadinstitute'` |
| Required APIs | `aiplatform.googleapis.com`, `iamcredentials.googleapis.com` |

**Service accounts (one per use case):**

| SA email | Use case | Workflow |
|---|---|---|
| `viral-ngs-cve-triage@viral-seq-ai.iam.gserviceaccount.com` | Weekly CVE triage | `.github/workflows/container-scan.yml` |

The pool + provider are reusable across use cases. Add a new SA per use case
(principle of least privilege; easier to audit and disable individually).

## Adding a New Claude-in-CI Use Case

You need GCP IAM Admin on `viral-seq-ai` for steps 1–3.

### 1. Create a service account scoped to the use case

```bash
gcloud iam service-accounts create <use-case-name> \
  --project=viral-seq-ai \
  --display-name="<Human-readable name>" \
  --description="Used by <workflow-file>.yml to invoke Claude on Vertex AI for <purpose>"
```

### 2. Grant minimum roles

Vertex invocation requires both:

```bash
gcloud projects add-iam-policy-binding viral-seq-ai \
  --member="serviceAccount:<use-case-name>@viral-seq-ai.iam.gserviceaccount.com" \
  --role="roles/aiplatform.user" \
  --condition=None

gcloud projects add-iam-policy-binding viral-seq-ai \
  --member="serviceAccount:<use-case-name>@viral-seq-ai.iam.gserviceaccount.com" \
  --role="roles/serviceusage.serviceUsageConsumer" \
  --condition=None
```

### 3. Bind the GitHub repo to the SA via WIF

```bash
PROJECT_NUMBER=$(gcloud projects describe viral-seq-ai --format='value(projectNumber)')

gcloud iam service-accounts add-iam-policy-binding \
  <use-case-name>@viral-seq-ai.iam.gserviceaccount.com \
  --project=viral-seq-ai \
  --role="roles/iam.workloadIdentityUser" \
  --member="principalSet://iam.googleapis.com/projects/${PROJECT_NUMBER}/locations/global/workloadIdentityPools/github-actions-pool/attribute.repository/broadinstitute/viral-ngs"
```

For tighter scope, use `attribute.ref` (branch) or `attribute.workflow` (specific
workflow file) in place of (or in addition to) `attribute.repository`.

### 4. Set GitHub repo variables

If the new use case lives alongside the CVE triage one, you can reuse
`GCP_PROJECT_ID` and `GCP_WIP_PROVIDER`. The SA email differs per use case —
either create a use-case-specific variable (e.g., `GCP_PR_REVIEW_SA_EMAIL`) or
hard-code it in the workflow file.

```bash
gh variable set GCP_<USECASE>_SA_EMAIL \
  --body "<use-case-name>@viral-seq-ai.iam.gserviceaccount.com" \
  --repo broadinstitute/viral-ngs
```

### 5. Add the workflow steps

See the canonical pattern below.

## Canonical Workflow Pattern

```yaml
permissions:
  id-token: write   # for OIDC token to GCP via WIF
  # plus whatever else your workflow needs (contents: read, issues: write, etc.)

steps:
  - name: Checkout repository
    uses: actions/checkout@v4
    with:
      # If Claude needs git history (`git log --grep`, `git show`), use 0.
      # Otherwise the default shallow clone is fine.
      fetch-depth: 0

  - name: Authenticate to GCP via Workload Identity Federation
    uses: google-github-actions/auth@v2
    with:
      workload_identity_provider: ${{ vars.GCP_WIP_PROVIDER }}
      service_account: ${{ vars.GCP_<USECASE>_SA_EMAIL }}

  - name: Claude on Vertex AI
    # Pin to commit SHA for supply-chain safety; bump when picking up new releases.
    uses: anthropics/claude-code-action@<FULL_40_CHAR_SHA>  # v1
    env:
      CLAUDE_CODE_USE_VERTEX: '1'
      CLOUD_ML_REGION: global              # Sonnet 4.6 supports the global endpoint
      ANTHROPIC_VERTEX_PROJECT_ID: ${{ vars.GCP_PROJECT_ID }}
    with:
      use_vertex: 'true'
      github_token: ${{ secrets.GITHUB_TOKEN }}    # avoids needing the Claude Code GitHub App
      claude_args: '--model claude-sonnet-4-6 --max-turns 30'
      settings: |
        {
          "permissions": {
            "allow": [
              "Read",
              "Write",
              "Bash(git log:*)",
              "Bash(grep:*)",
              "Bash(jq:*)"
              # ... add only what the prompt actually needs
            ]
          }
        }
      prompt: |
        <Your prompt here. Claude sees the GH workspace as cwd.>
```

## Gotchas (Things We Learned the Hard Way)

- **`iamcredentials.googleapis.com` must be enabled.** WIF impersonation calls
  this API; without it you get `IAM Service Account Credentials API has not
  been used in project ... before or it is disabled` from the Claude action,
  AFTER auth appears to succeed. Confusing.

- **Use `@v1`, not `@beta`.** `claude-code-action@beta` is the older API shape
  with separate `direct_prompt`, `max_turns`, `allowed_tools` inputs. `@v1`
  consolidates everything into `prompt` + `claude_args` + `settings`. Pin to a
  commit SHA, not the floating tag.

- **Both `use_vertex: 'true'` AND `CLAUDE_CODE_USE_VERTEX=1` are needed** —
  the input goes to the action wrapper, the env var goes to the underlying
  `claude-code` CLI.

- **Pass `github_token: ${{ secrets.GITHUB_TOKEN }}`** to skip the Claude
  Code GitHub App requirement. Without this you get
  `Error: Claude Code is not installed on this repository`.

- **The permission DSL doesn't take path globs on `Write` or paths on `Bash`
  command names.** `Write(/tmp/issues/**)` and `Bash(mkdir:/tmp/issues*)` are
  silently rejected and Claude hits `permission_denials_count > 0`. The
  `Bash(<cmd>:*)` pattern is for the *args after the command*, not paths
  embedded in the command name. For now, use unrestricted `Write` and rely
  on prompt instructions to constrain output paths.

- **Region `global` (recommended) gives dynamic routing across regions for
  Sonnet 4.6.** Pin to a specific region (e.g., `us-east5`, `europe-west1`)
  only if you need data-residency control.

- **`fetch-depth: 0` if Claude needs git history.** Default `actions/checkout`
  is shallow (depth=1); `git log --all --grep` and `git show <sha>` will
  produce empty results without unshallowing.

- **Always `jq`-query authoritative data sources in the prompt.** When Claude
  has multiple ways to learn a fact (training data vs reading a workspace
  file), tell it explicitly which is canonical and require the tool call.
  Otherwise it sometimes infers from training data.

## Cost / Safety

- **Cost gate:** invoke Claude only when the workflow has real work to do
  (e.g., new CVEs detected, PR opened by non-bot). Don't invoke on every
  schedule tick unconditionally.
- **Turn cap:** `--max-turns 30` is generous for triage-style tasks; tighten
  if you can. Observed: 6–10 turns for one-CVE analyses.
- **Cost order of magnitude:** Sonnet 4.6 ≈ $0.10–1 per non-trivial task
  (one-CVE analysis with full repo reading). Opus 4.7 is ~5× more expensive
  for marginal quality gain on most CI tasks.
- **Provider gate:** the OIDC provider attribute condition limits token
  minting to repos owned by `broadinstitute`. Other GitHub orgs cannot use
  this pool.
- **SA gate:** each SA's `workloadIdentityUser` binding limits which repos
  can impersonate it. Default scope: `attribute.repository/broadinstitute/viral-ngs`.
  Tighten with `attribute.ref` or `attribute.workflow` if needed.
- **Tool allowlist:** only allow tools the prompt actually uses. Avoid wildcard
  `Bash(*:*)` — name specific commands like `Bash(git log:*)`.

## Key Files

| File | Purpose |
|------|---------|
| `.github/workflows/container-scan.yml` | First use case (CVE triage); reference for the workflow pattern |
| `.agents/skills/container-vulns/SKILL.md` | The CVE triage playbook this infra serves |

## References

- [`anthropics/claude-code-action`](https://github.com/anthropics/claude-code-action)
- [Claude Code on Google Vertex AI](https://code.claude.com/docs/en/google-vertex-ai)
- [`google-github-actions/auth`](https://github.com/google-github-actions/auth) (WIF action)
