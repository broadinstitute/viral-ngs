#!/bin/sh
set -e

nosetests -v --with-xunit --with-coverage --nocapture \
    --cover-erase --cover-inclusive --cover-branches --cover-tests \
    --cover-package broad_utils,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util
