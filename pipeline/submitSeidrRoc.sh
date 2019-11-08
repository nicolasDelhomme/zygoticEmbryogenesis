#!/bin/bash
set -ex
# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
proj=u2019016
#indir=../data/seidr/results/backbone
indir=../data/seidr/results/aggregate
out=../data/seidr/roc/
gs=/mnt/picea/projects/functional-genomics-course/summer2018/goldStandard/Picea-abies_KEGG-based-positive-gold-standard.tsv

if [ ! -d $out ]; then
  mkdir -p $out
fi

# find the network files
#for f in $(find $indir -name "*.sf"); do
for f in $(find $indir -name "*.sf"); do
  fnam=$(basename ${f/.sf/})

  # identify networks
  read -r -a algos <<< $(seidr view -H $f | grep "\[A\]" | cut -d" " -f3 | xargs)
  
  # for every algos
  for i in $(seq 0 $(expr ${#algos[@]} - 1)); do
     echo $i
     sbatch -A $proj -o $out/${fnam}_roc.out \
     -e $out/${fnam}_roc.err \
     ../UPSCb-common/pipeline/runSeidrRoc.sh $f $gs $(expr $i + 1) \
     $out/${fnam}_${algos[$i]}_roc.tsv
  done
done
