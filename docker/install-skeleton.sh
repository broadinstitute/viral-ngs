#!/bin/bash

set -e -o pipefail

mkdir -p /opt/viral-ngs/viral-ngs-etc
ln -s /opt/viral-ngs/source /opt/viral-ngs/viral-ngs-etc/viral-ngs
ln /opt/docker/easy-deploy-viral-ngs.sh /opt/viral-ngs
ln /opt/docker/env_wrapper.sh /opt/viral-ngs

/opt/viral-ngs/easy-deploy-viral-ngs.sh setup-git-local
