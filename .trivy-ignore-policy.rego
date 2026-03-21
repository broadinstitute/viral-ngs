#
# bioinformatics-platform.rego
#
# Conservative Trivy ignore policy for bioinformatics container images
# running on genomics PaaS platforms as batch pipeline tasks.
#
# ASSUMPTIONS (document and verify these match your platform):
#   1. Containers run as non-interactive batch jobs (no shell sessions,
#      no web UIs, no Jupyter notebooks served from pipeline containers)
#   2. Containers have no inbound network listeners (no ports exposed)
#   3. Containers run with dropped capabilities (no CAP_SYS_ADMIN, etc.)
#   4. Containers do not run in privileged mode
#   5. Pipeline inputs are data files (FASTQ, BAM, VCF, reference genomes)
#      which may be untrusted or malformed
#
# If any assumption does not hold for a given image or use case,
# DO NOT apply this policy to that image.
#
# USAGE:
#   trivy image --ignore-policy bioinformatics-platform.rego \
#     --severity CRITICAL,HIGH --ignore-unfixed <image>
#
# IMPORTANT: Before using this policy, run your image with:
#   trivy image --format json <image> > scan.json
# and inspect the CVSS field structure. The field paths below
# (input.CVSS, input.CweIDs, etc.) reflect the Trivy JSON output
# structure. If your Trivy version uses different paths, adjust
# accordingly. The CVSS vector string location has changed across
# Trivy versions (see https://github.com/aquasecurity/trivy/issues/1627).
#
# CVSS VERSION SUPPORT:
#   This policy supports both CVSS v3.1 and CVSS v4.0 vector strings.
#   Trivy is transitioning to v4.0 for newer advisories. Some CVEs may
#   have only a v4.0 vector (no v3.1). The helper functions at the bottom
#   extract vectors from both versions and the rules are written to match
#   either format.
#
# VERSION: 2.0
# LAST REVIEWED: 2026-03-20
# REVIEW CADENCE: Quarterly, or when platform architecture changes
#

package trivy

default ignore = false

###############################################################################
# SECTION 1: PHYSICAL ACCESS REQUIRED (AV:P)
#
# Rationale: Cloud-hosted containers are never physically accessible.
# These CVEs require hands-on hardware interaction (USB, Firewire,
# JTAG, etc.) which is impossible in any cloud PaaS context.
#
# CVSS v3.1: AV:P
# CVSS v4.0: AV:P (same field name)
#
# Risk of false negative: Essentially zero. There is no scenario
# in which a pipeline container is physically accessible to an attacker.
# Confidence: Very High
###############################################################################

ignore {
    has_v3_field(input, "AV:P")
}

ignore {
    has_v4_field(input, "AV:P")
}

###############################################################################
# SECTION 2: ADJACENT NETWORK REQUIRED (AV:A)
#
# Rationale: Adjacent-network attacks require the attacker to be on
# the same physical or logical network segment (e.g., same VLAN,
# Bluetooth, local WiFi). In a cloud PaaS, pipeline containers run
# on orchestrated infrastructure where the attacker cannot place
# themselves on an adjacent segment.
#
# CVSS v3.1: AV:A
# CVSS v4.0: AV:A (same field name)
#
# Risk of false negative: Very low. Cloud networking abstractions
# make adjacent-network attacks impractical against pipeline containers.
# Confidence: Very High
###############################################################################

ignore {
    has_v3_field(input, "AV:A")
}

ignore {
    has_v4_field(input, "AV:A")
}

###############################################################################
# SECTION 3: USER INTERACTION REQUIRED (UI:R) + LOCAL VECTOR (AV:L)
#
# Rationale: Batch pipeline containers have no interactive user sessions.
# No human is clicking links, opening files in a GUI, or interacting
# with the container during execution. CVEs that require BOTH local
# access AND user interaction (e.g., tricking a user into opening a
# malicious file in a desktop app) are not exploitable in this context.
#
# NOTE: We require BOTH conditions (AV:L AND UI:R), not either alone.
# - AV:L alone is NOT safe to ignore (local privilege escalation
#   could be triggered by pipeline code)
# - UI:R alone is NOT safe to ignore for AV:N vulns (some network
#   vulns with UI:R involve clicking a link, which doesn't apply,
#   but others are ambiguous)
# - AV:L + UI:R together means "must have local access AND a human
#   must do something" - genuinely inapplicable in batch containers.
#
# CVSS v3.1: AV:L + UI:R
# CVSS v4.0: AV:L + UI:P (Passive) or UI:A (Active)
#   v4.0 splits "Required" into Passive (viewing content) and Active
#   (clicking/interacting). Both require a human, so both are safe to
#   ignore in batch containers.
#
# Risk of false negative: Very low for true batch pipeline containers.
# Confidence: High
###############################################################################

# v3: AV:L + UI:R
ignore {
    has_v3_field(input, "AV:L")
    has_v3_field(input, "UI:R")
}

# v4: AV:L + UI:P (Passive user interaction)
ignore {
    has_v4_field(input, "AV:L")
    has_v4_field(input, "UI:P")
}

# v4: AV:L + UI:A (Active user interaction)
ignore {
    has_v4_field(input, "AV:L")
    has_v4_field(input, "UI:A")
}

###############################################################################
# SECTION 4: HIGH PRIVILEGES REQUIRED + LOCAL VECTOR
#
# Rationale: CVEs requiring both local access and high (administrative)
# privileges assume the attacker already has elevated access to the
# system. In a properly configured container (non-root user, dropped
# capabilities), the process inside the container does not have high
# privileges to begin with. Combined with the local access requirement,
# this class of CVE is not practically exploitable.
#
# NOTE: We only ignore AV:L + PR:H, not AV:N + PR:H. A network-
# accessible vulnerability requiring high privileges may still be
# relevant if the service runs as a privileged user.
#
# CVSS v3.1: AV:L + PR:H
# CVSS v4.0: AV:L + PR:H (same field names)
#
# Risk of false negative: Low, assuming containers run as non-root.
# If your containers run as root, REMOVE THIS RULE.
# Confidence: High (conditional on non-root execution)
###############################################################################

ignore {
    has_v3_field(input, "AV:L")
    has_v3_field(input, "PR:H")
}

ignore {
    has_v4_field(input, "AV:L")
    has_v4_field(input, "PR:H")
}

###############################################################################
# SECTION 5: LOCAL ATTACK, SCOPE UNCHANGED (AV:L + S:U)
#
# Rationale: AV:L means the attacker already has local code execution
# inside the container. S:U (Scope Unchanged) means the impact does not
# cross a security boundary — it stays within the container.
#
# In an ephemeral batch container, this combination means: the attacker
# can already execute arbitrary code, and the vulnerability only lets
# them affect things inside a container that will be destroyed when the
# job completes. The vulnerability grants no capability the attacker
# does not already have.
#
# Contrast with S:C (Scope Changed): a local vulnerability that crosses
# the container-host boundary (e.g., container escape via kernel exploit)
# IS dangerous and is NOT ignored by this rule.
#
# CVSS v3.1: AV:L + S:U
# CVSS v4.0: AV:L + SC:N + SI:N + SA:N
#   v4.0 replaced the binary S:U/S:C with three subsequent-component
#   impact fields. SC:N + SI:N + SA:N means no impact on any component
#   beyond the vulnerable one — equivalent to v3's S:U.
#
# Risk of false negative: Low. The theoretical concern is that AV:L+S:U
# could include reading mounted secrets, but an attacker with code
# execution can already read those secrets directly.
# Confidence: High
###############################################################################

# v3: AV:L + S:U
ignore {
    has_v3_field(input, "AV:L")
    has_v3_field(input, "S:U")
}

# v4: AV:L + no subsequent-component impact
ignore {
    has_v4_field(input, "AV:L")
    has_v4_field(input, "SC:N")
    has_v4_field(input, "SI:N")
    has_v4_field(input, "SA:N")
}

###############################################################################
# SECTION 6: AVAILABILITY-ONLY IMPACT, SCOPE UNCHANGED
#
# Rationale: CVEs where the only impact is availability (DoS/resource
# exhaustion) and scope is unchanged mean: processing crafted input can
# crash or hang the affected process, but cannot read data (C:N),
# modify data (I:N), or affect other components (S:U).
#
# In ephemeral batch containers, a DoS means a single pipeline job
# fails or hangs until it hits its timeout or memory limit. This is
# operationally equivalent to a corrupted input file or OOM — the job
# fails, the container is destroyed, and the next job runs on a fresh
# container. There is no persistent state corruption, no data
# exfiltration, and no lateral movement.
#
# This rule applies regardless of attack vector (including AV:N),
# because the impact is strictly contained: even if triggered by
# network-delivered data, the worst outcome is one failed job.
#
# NOTE: This does NOT ignore:
#   - DoS with S:C / SC≠N / SI≠N / SA≠N (scope changed — could
#     affect host or other containers)
#   - DoS combined with any confidentiality or integrity impact
#     (C≠N or I≠N), which could indicate data leaks or corruption
#     alongside the crash
#
# CVSS v3.1: C:N + I:N + S:U (with any A value)
# CVSS v4.0: VC:N + VI:N + SC:N + SI:N + SA:N (with any VA value)
#
# Risk of false negative: Low. The concern would be if a DoS could be
# weaponized into a resource exhaustion attack against the compute
# platform (e.g., repeatedly submitting jobs with crafted inputs to
# burn credits). This is a business logic concern mitigated by job
# submission controls and cost alerts, not by container hardening.
# Confidence: High
###############################################################################

# v3: C:N + I:N + S:U (availability-only, scope unchanged)
ignore {
    has_v3_field(input, "C:N")
    has_v3_field(input, "I:N")
    has_v3_field(input, "S:U")
}

# v4: VC:N + VI:N + no subsequent-component impact (availability-only)
ignore {
    has_v4_field(input, "VC:N")
    has_v4_field(input, "VI:N")
    has_v4_field(input, "SC:N")
    has_v4_field(input, "SI:N")
    has_v4_field(input, "SA:N")
}

###############################################################################
# HELPER FUNCTIONS: Extract and match CVSS vector strings
#
# Trivy's JSON structure nests CVSS data under input.CVSS with vendor
# keys. The vector string location varies by data source. We check
# multiple common paths and prefer NVD.
#
# CVSS v3.1 vectors look like: CVSS:3.1/AV:N/AC:L/PR:N/UI:N/S:U/C:N/I:N/A:H
# CVSS v4.0 vectors look like: CVSS:4.0/AV:N/AC:L/AT:N/PR:N/UI:N/VC:N/VI:N/VA:H/SC:N/SI:N/SA:N
#
# Fields are slash-delimited key:value pairs. The has_vX_field helpers
# check for a field anywhere in the vector, handling both mid-string
# (/field/) and end-of-string (/field) positions.
#
# IMPORTANT: Run `trivy image --format json <your-image>` and inspect
# the .Vulnerabilities[].CVSS structure to confirm these paths work
# for your Trivy version. If the structure differs, update these
# functions accordingly.
###############################################################################

# --- CVSS v3.1 vector extraction ---

get_v3_vector(vuln) = vector {
    vector := vuln.CVSS.nvd.V3Vector
} else = vector {
    vector := vuln.CVSS.redhat.V3Vector
} else = vector {
    vector := vuln.CVSS.ghsa.V3Vector
} else = vector {
    some vendor
    vector := vuln.CVSS[vendor].V3Vector
} else = "" {
    true
}

# --- CVSS v4.0 vector extraction ---

get_v4_vector(vuln) = vector {
    vector := vuln.CVSS.nvd.V40Vector
} else = vector {
    vector := vuln.CVSS.redhat.V40Vector
} else = vector {
    vector := vuln.CVSS.ghsa.V40Vector
} else = vector {
    some vendor
    vector := vuln.CVSS[vendor].V40Vector
} else = "" {
    true
}

# --- Field matching helpers ---
# Check if a CVSS vector contains a specific field value.
# Handles both mid-string (/AV:N/) and end-of-string (/AV:N) positions.

has_v3_field(vuln, field) {
    cvss_vector := get_v3_vector(vuln)
    cvss_vector != ""
    contains(cvss_vector, concat("", ["/", field, "/"]))
}

has_v3_field(vuln, field) {
    cvss_vector := get_v3_vector(vuln)
    cvss_vector != ""
    endswith(cvss_vector, concat("", ["/", field]))
}

has_v4_field(vuln, field) {
    cvss_vector := get_v4_vector(vuln)
    cvss_vector != ""
    contains(cvss_vector, concat("", ["/", field, "/"]))
}

has_v4_field(vuln, field) {
    cvss_vector := get_v4_vector(vuln)
    cvss_vector != ""
    endswith(cvss_vector, concat("", ["/", field]))
}

###############################################################################
# RULES INTENTIONALLY NOT INCLUDED (and why):
#
# 1. AV:N (Network attack vector) — NOT blanket-ignored.
#    Even though batch pipeline containers typically have no inbound
#    listeners, some AV:N CVEs involve outbound connections triggered
#    by processing attacker-influenced data (e.g., Log4Shell). We
#    cannot safely blanket-ignore network-vector CVEs. However,
#    Section 6 does ignore AV:N CVEs that are availability-only with
#    no scope change, since the worst outcome is a crashed job.
#
# 2. UI:R alone (without AV:L) — NOT ignored.
#    Some AV:N + UI:R vulnerabilities involve scenarios like processing
#    a crafted file that triggers a callback, which could be borderline
#    relevant if pipeline inputs are not fully trusted.
#
# 3. CWE-based filters — NOT included by default.
#    Filtering by CWE class (e.g., ignoring all CWE-79 XSS in a
#    non-web context) is tempting but risky as a default policy.
#    CWE classifications can be inaccurate, and a single misclassified
#    CVE could slip through. If your platform has specific architectural
#    mitigations (e.g., provably no web server in any image), you may
#    add CWE rules, but this should be a deliberate per-platform choice.
#
# 4. Specific CVE IDs — NOT included.
#    Use a .trivyignore.yaml or VEX document for individual CVE
#    exceptions with per-CVE justifications. The Rego policy should
#    capture architectural/class-level mitigations, not one-off
#    exceptions.
#
# 5. AV:L + S:C — NOT ignored.
#    Local vulnerabilities with Scope Changed (S:C) can cross security
#    boundaries (e.g., container escape via kernel exploit). These are
#    dangerous even in ephemeral containers. Only AV:L + S:U (Scope
#    Unchanged) is ignored — see Section 5 above.
#
# 6. Inbound-listener-only server CVEs — NOT categorically ignored.
#    Many AV:N CVEs in fat JARs (Jetty, ZooKeeper, Netty server-side)
#    require an active network listener that we never start. However,
#    CVSS does not distinguish inbound-listener vs. data-processing
#    attack surfaces within AV:N. Adding package-name-based exceptions
#    here would be fragile and is better handled in .trivyignore with
#    per-CVE justification documenting that the server component is
#    never instantiated.
###############################################################################
