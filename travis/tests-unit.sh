#!/bin/bash
set -e

echo travis_fold:start:tests-unit

nosetests -v --with-xunit --with-coverage --nocapture \
    --cover-inclusive --cover-branches --cover-tests \
    --cover-package broad_utils,illumina,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util \
    -w test/unit/

echo travis_fold:end:tests-unit
