## Differential gene expression scripts

Scripts that I wrote to assess the extent of differential gene expression in 
*Populus balsamifera*, and to estimate rates of molecular evolution on 
differentially expressed genes. 
 
### Contents

1. [README.md](https://github.com/BrianSanderson/gene-expression/blob/master/R/README.md): this file

2. [runDGE.R](https://github.com/BrianSanderson/gene-expression/blob/master/R/runDGE.R) An R function that will estimate differential gene expression between two groups using the consensus of three different metrics: DeSeq2, limma, and edgeR.

3. [getSpecificGenes.R](https://github.com/BrianSanderson/gene-expression/blob/master/R/getSpecificGenes.R) An R function that will quantify sex- and tissue-specific patterns of gene expression. This code really is a slog. Please feel free to reach out if you have a more elegant solution.

4. [bestBlastHits.R](https://github.com/BrianSanderson/gene-expression/blob/master/R/bestBlastHits.R) An R script that I wrote to select the best homolog for *P. trichocarpa* sequences in *P. tremula* and *P. euphratica* based on BLAST results. See comments in the header to run from the command line using Rscript 
