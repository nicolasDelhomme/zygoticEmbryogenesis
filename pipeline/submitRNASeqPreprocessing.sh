#!/bin/bash -l

set -eux

proj=u2019016
mail=mist0069@student.umu.se

in=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/raw/RNA-Seq
out=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/RNA-Seq
start=2
end=6

##corrected module calls to match
module load bioinfo-tools FastQC Trimmomatic sortmerna

for f in $(find $in -name "*_1.fastq.gz"); 
do
  fnam=$(basename ${f/_1.fastq.gz/})
  bash ../UPSCb-common/pipeline/runRNASeqPreprocessing.sh -s $start -e $end \
  $proj $mail $f $in/${fnam}_2.fastq.gz $out
done

