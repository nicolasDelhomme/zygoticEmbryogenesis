# A transcriptome atlas of zygotic and somatic embryogenesis in Norway spruce

This repository is about the spruce zygotic embryogenesis analysis

## Authors
Katja Stojkovič, Camilla Canovi, Kim-Cuong Le, Iftikhar Ahmad, Ioana Gaboreanu, Sofie Johansson, Nicolas Delhomme, Ulrika Egertsdotter, Nathaniel R. Street

## Abstract

Somatic embryogenesis (SE) is a powerful model system for studying embryo development and an important method for scaling up availability of elite and climate-adapted genetic material of Norway spruce (Picea abies L. Karst). However, there are several steps during the development of the somatic embryo (Sem) that are suboptimal compared to zygotic embryo (Zem) development. These differences are poorly understood and result in substantial yield losses during plant production, which limits cost-effective large-scale production of SE plants. This study presents a comprehensive data resource profiling gene expression during zygotic and somatic embryo development to support studies aiming to advance understanding of gene regulatory programmes controlling embryo development. Transcriptome expression patterns were analysed during zygotic embryogenesis (ZE) in Norway spruce, including separated samples of the female gametophytes and Zem, and at multiple stages during SE. Expression data from eight developmental stages of SE, starting with pro-embryogenic masses (PEMs) up until germination, revealed extensive modulation of the transcriptome between the early and mid-stage maturing embryos and at the transition of desiccated embryos to germination. Comparative analysis of gene expression changes during ZE and SE identified differences in the pattern of gene expression changes and functional enrichment of these provided insight into the associated biological processes. Orthologs of transcription factors known to regulate embryo development in angiosperms were differentially regulated during Zem and Sem development and in the different zygotic embryo tissues, providing clues to the differences in development observed between Zem and Sem. This resource represents the most comprehensive dataset available for exploring embryo development in conifers.

## Links

[Manuscript](https://doi.org/10.1111/tpj.17087])

## Citation

This repository: [![DOI](https://zenodo.org/badge/940495767.svg)](https://doi.org/10.5281/zenodo.14944673)

## Setup
```{bash setup, eval=FALSE}
ln -s /mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series data
ln -s /mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/analysis analysis
```

## Differential expression of genes in zygotic embryogenesis (Katja)
Stages of zygotic embryogenesis were renamed from B[0-9]+ to Z[0-9]+ to indicate they come from zygotic embryogenesis, compared to stages named S[0-9], coming from somatic embryogenesis, to which they were compared.

Biological quality assessment:
exclude B10, as it has only one replicate
src/R/BiologicalQC-ZE-FMG.R

DE:
src/R/DifferentialExpression-ZE-FMG-allStages.R
results in analysis/ZE-FMG-allStages_dupl_Samples
(compare and) delete analysis/ZE_FMG_designReplicate?

Exploration of DEGs:
src/R/ExploringDEGs-ZE-FMG-SE.R

