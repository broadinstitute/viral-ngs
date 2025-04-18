#!/bin/bash
set -e -o pipefail

echo "which gatk3: $(which gatk3)"
echo "gatk3 version: $(gatk3 --version)"

# if the gatk3 command is available
if [[ ! -z "$(which gatk3)" ]]; then
    echo "'gatk3' command found, symlinking to 'gatk'..."
    ln -sf "$(which gatk3)" "$(dirname $(which gatk3))/gatk"
else
    echo "'gatk3' command not found, skipping symlink creation."
fi


