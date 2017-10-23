#!/bin/bash
set -e -x -o pipefail

# This script performs various packing and deployment operations.
# It assumes it will be caused as a deploy hook of TravisCI; ex.:
#
# deploy:
#   provider: script
#   script: travis/deploy.sh $TRAVIS_TAG
#   on:
#     tags: true
#     all_branches: master


# way to get the absolute path to this script that should
# work regardless of whether or not this script has been sourced
# Find original directory of bash script, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        SOURCE="$(readlink "$SOURCE")"
    else
        SOURCE="$(readlink -f "$SOURCE")"
    fi
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT=$SOURCE
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$SCRIPT_DIRNAME" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"

PKG_VERSION=$1
CONDA_PACKAGE_OUTDIR=packacking/conda-packages

# === Build conda package

echo "Python binary: $(which python)"
echo "Python version: $(python --version)"

conda config --set anaconda_upload yes

if [ $BUILD_PACKAGE = "true" ]; then
    # If this is a PR, on the master branch, or is a tag, render and build the conda package. If it is a tag, also upload to anaconda.org
    if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then
        echo "Rendering and building conda package..."
        # Render recipe from template and dependency files, setting the tag as the current version
        # if this is a tag build+upload, otherwise just test building
        if [ -n "$TRAVIS_TAG" ]; then
             # if the ANACONDA_TOKEN is defined (not on an external branch)
            if [ ! -z "$ANACONDA_TOKEN" ]; then
                
                python packaging/conda-recipe/render-recipe.py "$PKG_VERSION" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt && \
                CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" packaging/conda-recipe/viral-ngs && \
                
                REPO=broadinstitute/viral-ngs
                TAG=$(if [ "$TRAVIS_TAG" == "master" ]; then echo "latest"; else echo "$TRAVIS_TAG" ; fi)
                VIRAL_NGS_VERSION=$(echo "$TRAVIS_TAG" | perl -lape 's/^v(.*)/$1/g') # strip 'v' prefix
                tar -czh . | docker build --build-arg VIRAL_NGS_VERSION=$VIRAL_NGS_VERSION --rm -t "$REPO:$VIRAL_NGS_VERSION" - | tee >(grep "Successfully built" | perl -lape 's/^Successfully built ([a-f0-9]{12})$/$1/g' > build_id) | grep ".*" && build_image=$(head -n 1 build_id) && rm build_id

                if [[ ! -z $build_image ]]; then
                    echo "build_image: $build_image"

                    # export and import to squash an image to one layer. 
                    # only containers can be exported, so we need to make a temp one
                    #temp_container_id=$(docker create "$REPO:$VIRAL_NGS_VERSION-build")
                    #docker export "$temp_container_id" | docker import - "$REPO:$VIRAL_NGS_VERSION-run-precursor"
                    #echo "FROM $REPO:$VIRAL_NGS_VERSION-run-precursor"'
                    #      ENTRYPOINT ["/opt/viral-ngs/env_wrapper.sh"]' | docker build -t "$REPO:$VIRAL_NGS_VERSION" - | tee >(grep "Successfully built" | perl -lape 's/^Successfully built ([a-f0-9]{12})$/$1/g' > build_id) | grep ".*" && build_image=$(head -n 1 build_id) && rm build_id

                    docker run --rm $build_image illumina.py && docker login -u "$DOCKER_USER" -p "$DOCKER_PASS" && docker push "$REPO:$VIRAL_NGS_VERSION"
                else
                    echo "Docker build failed."
                    exit 1
                fi
            else
                echo "ANACONDA_TOKEN is not defined. Conda package upload is only supported for branches on the original repository."
            fi
        else
            # make a directory to hold the built conda package
            mkdir -p CONDA_PACKAGE_OUTDIR
            # build the conda package
            described_version="$(./assembly.py --version)"
            package_suffix="+$described_version"
            python packaging/conda-recipe/render-recipe.py "0.0.0$package_suffix" --package-name "viral-ngs-dev" --download-filename "$TRAVIS_BRANCH" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt && \
            #CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --no-anaconda-upload --output-folder "$CONDA_PACKAGE_OUTDIR" packaging/conda-recipe/viral-ngs
            CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --output-folder "$CONDA_PACKAGE_OUTDIR" packaging/conda-recipe/viral-ngs
            cp $CONDA_PACKAGE_OUTDIR/*/*.tar.bz2 ./docker/
            # build the docker image, and try to run it
            rc="$?"
            if [[ "$rc" == "0" ]]; then
                pushd ./docker

                REPO=broadinstitute/viral-ngs-dev
                TAG=$described_version

                tar -czh . | docker build --rm - | tee >(grep "Successfully built" | perl -lape 's/^Successfully built ([a-f0-9]{12})$/$1/g' > build_id) | grep ".*" && build_image=$(head -n 1 build_id) && rm build_id
                popd
                if [[ ! -z $build_image ]]; then
                    echo "build_image: $build_image"
                    docker run --rm $build_image illumina.py || exit $?
                else
                    echo "Docker build failed."
                    exit 1
                fi
            fi
        fi
    fi  
else
    echo "Not building a package for this slot in the build matrix."
fi
