#!/bin/bash
set -e -o pipefail

cached_fetch_jar_from_github () {
	_github_org=$1
	_repo_name=$2
	_tool_name=$3
	_jar_version=$4
	_jar_fname="$_tool_name-$_jar_version.jar"
	if [ ! -f $CACHE_DIR/$_jar_fname ]; then
		echo "Fetching $_jar_fname"
		wget --quiet https://github.com/$_github_org/$_repo_name/releases/download/$_jar_version/$_jar_fname
		mv $_jar_fname $CACHE_DIR
	else
		echo "Using cached $_jar_fname"
	fi
	ln -s $CACHE_DIR/$_jar_fname $_tool_name.jar
}

cached_fetch_jar_from_github broadinstitute cromwell womtool 33.1
cached_fetch_jar_from_github broadinstitute cromwell cromwell 33.1
cached_fetch_jar_from_github dnanexus dxWDL dxWDL 0.72

TGZ=dx-toolkit-v0.285.0-ubuntu-16.04-amd64.tar.gz
if [ ! -f $CACHE_DIR/$TGZ ]; then
	echo "Fetching $TGZ"
	wget --quiet https://dnanexus-sdk.s3.amazonaws.com/$TGZ
	mv $TGZ $CACHE_DIR
else
	echo "Using cached $TGZ"
fi
tar -xzpf $CACHE_DIR/$TGZ

#echo "Fetching quayctl"
#wget --quiet https://github.com/coreos/quayctl/releases/download/v0.0.1/quayctl-linux-x64
#mv quayctl-linux-x64 quayctl
#chmod a+x quayctl
