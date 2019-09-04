#!/bin/bash

#set -ex

proj=b2013162
mail=izabela.dobrowolska@slu.se
in="/proj/$proj/nobackup/sRNA-SE-germinants/fastqc/reaper"
out="/proj/$proj/nobackup/sRNA-SE-germinants/multiqc/reaper"

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools MultiQC

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err -A $proj $UPSCb/pipeline/runMultiQC.sh $in $out
