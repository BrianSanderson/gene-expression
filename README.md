## Differential gene expression scripts

Scripts that I wrote to assess the extent of differential gene expression in *Populus balsamifera*, and to estimate rates of molecular evolution on differentially expressed genes.
 
### Contents

1. [README.md](https://github.com/BrianSanderson/genomics-scripts/blob/master/gene-expression/README.md): this file

2. [subsetBlastResults.R](https://github.com/BrianSanderson/genomics-scripts/blob/master/gene-expression/subsetBlastHits.R) An R script that I wrote to try to select the best homolog for *P. trichocarpa* sequences in *P. tremula* and *P. euphratica* based on BLAST results. It can be run from the command line using Rscript, with the first argument referencing BLAST output with the following fields: "qseqid", "qlen", "sseqid", "length", "pident", "evalue", "bitscore."
