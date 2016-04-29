#!/bin/bash
#set -e

nosetests -v \
    --logging-clear-handlers \
    --with-timer \
    --with-xunit \
    -w test/unit/
    # --with-xunit --with-coverage \
    # --cover-inclusive --cover-branches --cover-tests \
    # --cover-package broad_utils,illumina,assembly,interhost,intrahost,metagenomics,ncbi,read_utils,reports,taxon_filter,tools,util \
