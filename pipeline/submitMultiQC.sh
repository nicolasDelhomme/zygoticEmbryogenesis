#!/bin/bash

#set -ex

proj=u2019016
mail=mist0069@student.umu.se
in=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/RNA-Seq
out=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/RNA-Seq/multiqc

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools multiqc

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err -A $proj ../UPSCb-common/pipeline/runMultiQC.sh $in $out
