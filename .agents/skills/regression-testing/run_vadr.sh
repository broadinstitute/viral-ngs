#!/bin/bash
# Run VADR on a single FASTA file.
# Inputs (env vars set by dsub):
#   FASTA      - input FASTA file (localized by dsub)
#   VADR_OPTS  - vadr options string
#   MIN_LEN    - minimum sequence length (optional)
#   MAX_LEN    - maximum sequence length (optional)
#   MODEL_URL  - URL to vadr model tarball
#   MODEL_SUB  - subdirectory within model tarball (optional)
# Outputs (env vars set by dsub):
#   NUM_ALERTS - file to write num_alerts integer
#   ALERTS_TSV - file to write alerts TSV
#   VADR_TGZ   - file to write full vadr output tarball

set -e

BASENAME=$(basename "${FASTA}" .fasta)

# Download and unpack VADR models
if [ -n "${MODEL_URL}" ]; then
  mkdir -p vadr-untar
  curl -fsSL "${MODEL_URL}" | tar -C vadr-untar -xzf -
  ln -s vadr-untar/*/ vadr-models
else
  ln -s /opt/vadr/vadr-models vadr-models
fi

if [ -n "${MODEL_SUB}" ]; then
  VADR_MODEL_DIR="vadr-models/${MODEL_SUB}"
else
  VADR_MODEL_DIR="vadr-models"
fi

# Build trim args
TRIM_ARGS=""
if [ -n "${MIN_LEN}" ]; then
  TRIM_ARGS="${TRIM_ARGS} --minlen ${MIN_LEN}"
fi
if [ -n "${MAX_LEN}" ]; then
  TRIM_ARGS="${TRIM_ARGS} --maxlen ${MAX_LEN}"
fi

# Remove terminal ambiguous nucleotides
/opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
  "${FASTA}" ${TRIM_ARGS} > "${BASENAME}.trimmed.fasta"

# Run VADR
v-annotate.pl \
  ${VADR_OPTS} \
  --split --cpu $(nproc) \
  --mdir "${VADR_MODEL_DIR}" \
  "${BASENAME}.trimmed.fasta" \
  "${BASENAME}"

# Package outputs
tar -C "${BASENAME}" -czf "${VADR_TGZ}" .

# Extract alerts
cat "${BASENAME}/${BASENAME}.vadr.alt.list" | cut -f 5 | tail -n +2 > alerts.tsv
cp alerts.tsv "${ALERTS_TSV}"
wc -l < alerts.tsv > "${NUM_ALERTS}"

echo "VADR complete. Alerts: $(cat ${NUM_ALERTS})"
