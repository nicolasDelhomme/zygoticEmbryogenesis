#!/bin/bash -l

## be verbose and print
set -eux

proj=u2019016
mail=mist0069@student.umu.se

## process the argument
in=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all.phase.gff3.CDS.fa
out=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/salmon/Pabies1.0-all-phase.gff3.CDSandLTR-TE


  
  ## execute
sbatch -A $proj --mail-user=$mail \
  -e $out/index.err -o $out/index.out -J salmonIndex \
  ../UPSCb-common/pipeline/runSalmonIndex.sh $in $out