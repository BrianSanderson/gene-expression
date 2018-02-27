## Differential gene expression scripts

Scripts that I wrote to assess the extent of differential gene expression in 
*Populus balsamifera*, and to estimate rates of molecular evolution on 
differentially expressed genes.
 
### Contents

1. [README.md](https://github.com/BrianSanderson/genomics-scripts/blob/master/gene-expression/README.md): this file

2. [DGEfunction.R](https://github.com/BrianSanderson/genomics-scripts/blob/master/gene-expression/DGEfunction.R) An R function that will estimate differential
gene expression between two groups using the consensus of three different
metrics: DeSeq2, limma, and edgeR.

3. [subsetBlastResults.R](https://github.com/BrianSanderson/genomics-scripts/blob/master/gene-expression/subsetBlastHits.R) An R script that I wrote to 
select the best homolog for *P. trichocarpa* sequences in *P. tremula* 
and *P. euphratica* based on BLAST results. See comments in the header to 
run from the command line using Rscript 
