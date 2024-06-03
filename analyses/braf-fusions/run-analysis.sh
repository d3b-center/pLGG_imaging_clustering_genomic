#!/bin/bash

set -e
set -o pipefail

# run breakpoint by imaging cluster
Rscript -e "rmarkdown::render('01-fusion-breakpoints-by-imaging-cluster.Rmd')"

# run survival by imaging cluster
Rscript -e "rmarkdown::render('02-imaging-cluster-survival.Rmd')"