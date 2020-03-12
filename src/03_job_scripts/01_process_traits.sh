#!/bin/bash

for script in src/03_job_scripts/01_helpers/*; do
  name=$(basename $script | sed 's/.R//')
  if [[ ! -d "analyses/processed_traits/$name" ]]; then
    Rscript $script
  fi
done
