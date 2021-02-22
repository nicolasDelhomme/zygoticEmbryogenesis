#!/bin/bash

script=~/Git/zygoticEmbryogenesis/src/R/ExploringDEGs-ZE-FMG-SE_runEST.R
dir=~/Git/zygoticEmbryogenesis/analysis/DE/ZE-FMG-allStages_duplSsamples

sbatch -A u2019016 -t 2-00:00:00 -o $dir/ExploringDEGs_runEST.out -e $dir/ExploringDEGs_runEST.err $UPSCb/UPSCb-common/pipeline/runRmarkdown.sh $script