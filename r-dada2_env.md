For running dada2 R on premise

Create R env

```
module load anaconda/colsa

conda create -n r-dada2_env r-essentials r-base
conda activate r-dada2_env
```

Because of issue referenced [here](https://github.com/benjjneb/dada2/issues/417) an additional comiler

```
R
>if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
>BiocManager::install("dada2")
```
### The above works as of 02/18/2020. Some earlier versions of these programs had bug at install. See below for potential fixes if running into problems...


Because of issue referenced [here](https://github.com/benjjneb/dada2/issues/417) an additional comiler

```
R
>if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
>BiocManager::install("dada2", version = "3.9")
```
This is currently failing on install of Rhtslib. It looks the install of this package should be fixed in a couple of days (from 9/18/19), but for now comments [here](https://github.com/Bioconductor/Rhtslib/issues/9)
```
wget https://bioconductor.org/packages/release/bioc/src/contrib/Rhtslib_1.16.1.tar.gz
tar xvzf Rhtslib_1.16.1.tar.gz
cd Rhtslib/src/htslib-1.7/
```
In Makefile and Makefile.rhtslib comment the lines `CPPFLAGS =` and `LDFLAGS =`
Change the `CFLAGS =` to `CFLAGS +=`

```
cd
tar -czvf Rhtslib_1.16.1.tar.gz ./Rhtslib/

R CMD INSTALL Rhtslib_1.16.1.tar.gz
```
then
```
R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")
```
This still fails on ShortRead package. try:
```
# Make sure channels are configured correctly
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Install GenomicAlignments
conda install bioconductor-genomicalignments
conda install bioconductor-shortread
```
then
```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")
```

