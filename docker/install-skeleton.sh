#!/bin/bash

set -e -o pipefail

mkdir -p /opt/viral-ngs/viral-ngs-etc
ln -s /opt/viral-ngs/source /opt/viral-ngs/viral-ngs-etc/viral-ngs
ln -s /opt/docker/easy-deploy-viral-ngs.sh /opt/viral-ngs
ln -s /opt/docker/env_wrapper.sh /opt/viral-ngs

/opt/docker/easy-deploy-viral-ngs.sh setup-git-local
