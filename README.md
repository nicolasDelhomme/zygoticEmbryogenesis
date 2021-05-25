# zygoticEmbryogenesis
Spruce zygotic embryogenesis 

## Setup
```{bash setup, eval=FALSE}
ln -s /mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series data
/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/analysis analysis
```

## Differential expression of genes in zygotic embryogenesis (Katja)
Stages of zygotic embryogenesis were renamed from B[0-9]+ to Z[0-9]+ to indicate they come from zygotic embryogenesis, compared to stages named S[0-9], coming from somatic embryogenesis, to which they were compared.

Biological quality assesment:
exclude B10, as it has only one replicate
src/R/BiologicalQC-ZE-FMG.R

DE:
src/R/DifferentialExpression-ZE-FMG-allStages.R
results in analysis/ZE-FMG-allStages_dupl_Samples
(compare and) delete analysis/ZE_FMG_designReplicate?

Exploration of DEGs:
src/R/ExploringDEGs-ZE-FMG-SE.R

