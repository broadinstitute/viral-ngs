#!/bin/bash

set -e

# upgrade docker (Travis's is old)
sudo apt-get update
sudo apt-get -y -o Dpkg::Options::="--force-confnew" install docker-ce
docker --version
