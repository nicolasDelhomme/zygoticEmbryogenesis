#!/bin/bash

# vars
account=u2019016
mail=nicolas.delhomme@slu.se
in=../data/seidr/results/sf
out=../data/seidr/results/aggregate

# modules
module load bioinfo-tools seidr-devel

# dir
if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
sbatch -A $account --mem=128GB --mail-user=$mail \
-e $out/aggregate.err -o $out/aggregate.err \
-J ZE-aggregate ../UPSCb-common/pipeline/runSeidrAggregate.sh $out \
$(find $in -name "*.sf" -exec realpath "{}" \;)

