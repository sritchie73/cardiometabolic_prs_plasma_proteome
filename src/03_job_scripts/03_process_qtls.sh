#!/bin/bash

for script in src/03_job_scripts/03_helpers/*; do
  Rscript $script
done
