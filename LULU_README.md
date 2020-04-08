LULU ASV filtering (Froslev et al. 2017, Nature Communications, 8(1), 1188). This is performed on files that were processed for read quality and ITS extracted. Files from separarte runs were then concatenated by sample and DADA2 ASV calling with `pool = T` was performed.

for install `https://github.com/tobiasgf/lulu`

First perform sequence similarity matching

```
conda activate itsxpress_3p
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
mkdir LULU

vsearch --usearch_global dada2_out/ASVs.fa --db dada2_out/ASVs.fa --self --id .84 --iddef 1 --userout LULU/match_list.0.84min.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

conda deactivate
```

Run LULU algorithm at min_similarity 84, 90, 93, 95
```
Rscript ~/repo/neonectria_barcoding_012220/LULU/LULU.r
```
Count the number of ASVs with species and genus names reassigned. Specieas are counted as reassigned if the child is not NA (but is named) and is assigned to something other than the named species (including if it is then assigned to NA). Genus is counted the same, is also counted even if species level naming is already counted.
```
sims=( 84 90 93 95 )
for i in ${sims[@]}
do(
echo $i
perl ~/repo/neonectria_barcoding_012220/LULU/count_species_and_genera_reassigned.pl LULU/merged_asvs_${i}.txt LULU/which_reassigned_${i}.txt
)
done
```

