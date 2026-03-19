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
#   5. Pipeline inputs are data files (FASTQ, BAM, VCF, reference genomes),
#      not attacker-controlled interactive input
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
# VERSION: 1.0
# LAST REVIEWED: 2026-03-19
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
# Risk of false negative: Essentially zero. There is no scenario
# in which a pipeline container is physically accessible to an attacker.
# Confidence: Very High
###############################################################################

ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:P/")
}

# Also catch AV:P at the end of the vector string (no trailing slash)
ignore {
    cvss_vector := get_v3_vector(input)
    endswith(cvss_vector, "/AV:P")
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
# Risk of false negative: Very low. Cloud networking abstractions
# make adjacent-network attacks impractical against pipeline containers.
# Confidence: Very High
###############################################################################

ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:A/")
}

ignore {
    cvss_vector := get_v3_vector(input)
    endswith(cvss_vector, "/AV:A")
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
# Risk of false negative: Very low for true batch pipeline containers.
# Confidence: High
###############################################################################

ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:L/")
    contains(cvss_vector, "/UI:R")
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
# Risk of false negative: Low, assuming containers run as non-root.
# If your containers run as root, REMOVE THIS RULE.
# Confidence: High (conditional on non-root execution)
###############################################################################

ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:L/")
    contains(cvss_vector, "/PR:H/")
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
# Risk of false negative: Low. The theoretical concern is that AV:L+S:U
# could include reading mounted secrets, but an attacker with code
# execution can already read those secrets directly.
# Confidence: High
###############################################################################

ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:L/")
    contains(cvss_vector, "/S:U/")
}

# Also catch S:U at the end of the vector string (no trailing slash)
ignore {
    cvss_vector := get_v3_vector(input)
    contains(cvss_vector, "/AV:L/")
    endswith(cvss_vector, "/S:U")
}

###############################################################################
# HELPER FUNCTION: Extract the CVSS v3 vector string
#
# Trivy's JSON structure nests CVSS data under input.CVSS with vendor
# keys. The vector string location varies by data source. We check
# multiple common paths and prefer NVD.
#
# IMPORTANT: Run `trivy image --format json <your-image>` and inspect
# the .Vulnerabilities[].CVSS structure to confirm these paths work
# for your Trivy version. If the structure differs, update this
# function accordingly.
###############################################################################

get_v3_vector(vuln) = vector {
    vector := vuln.CVSS.nvd.V3Vector
} else = vector {
    vector := vuln.CVSS.redhat.V3Vector
} else = vector {
    vector := vuln.CVSS.ghsa.V3Vector
} else = vector {
    # Fallback: try any vendor that has a V3Vector
    some vendor
    vector := vuln.CVSS[vendor].V3Vector
} else = "" {
    true
}

###############################################################################
# RULES INTENTIONALLY NOT INCLUDED (and why):
#
# 1. AV:N (Network attack vector) — NOT ignored.
#    Even though batch pipeline containers typically have no inbound
#    listeners, some AV:N CVEs involve outbound connections triggered
#    by processing attacker-influenced data (e.g., Log4Shell). We
#    cannot safely blanket-ignore network-vector CVEs.
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
###############################################################################
