#!/bin/bash
set -e

nosetests -v --with-xunit --with-coverage \
    --cover-inclusive --cover-branches --cover-tests \
    --cover-package broad_utils,illumina,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util \
    -w test/unit/
