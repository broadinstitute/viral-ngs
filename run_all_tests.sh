#!/bin/sh

 set -e -x -o pipefail

#python -m unittest discover -v
# since test.unit.test_tools.TestToolsInstallation uses unittest-incompatible
# test generators, we have to test using nosetests
nosetests -v --with-xunit --with-coverage \
    --cover-inclusive --cover-branches --cover-tests \
    --cover-package broad_utils,illumina,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util \
    test/unit/ test/integration/
