For running dada2 R on premise

Create conda env. There are some additional flags here that are required to avoid Rcpp package version that result in sefaults on multithreading https://github.com/benjjneb/dada2/issues/684

```
module load anaconda/colsa

 conda create -n dada2-check -c conda-forge -c bioconda -c defaults --override-channels bioconductor-dada2
 ```
