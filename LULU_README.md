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
At 93% `minimum_similarity` there is minimal taxonomic reassignment at the genus and species level, but a similar percentage of ASVs are "discarded" (11.3%) as the default (84% seq sim) cutoff. Range of lumped/discarded ASVs is 9.9%-14% for 95-84% seq. sim., respectively. See table `LULU_discarded_ASVs.xlsx` for the numbers of taxa reassigned at genus and species level at different cutoffs. In addition, the one species reassigned to a different ASV at 93% cutoff that is also named at the species level at seems unlikely to be observed in North American beech trees. This ASV was named as *Fusarium algeriense* which has previously observed in wheat crops in Algeria, and is a relatively poor match to the identified species (UNITE-BLAST of the sequence does not turn up Fusarium algeriense as a potential match).

At 90% sequence similarity additional ASVs are assigned different names, including an ASV that has a nearly 100% sequence similarity match to a named species in UNITE, and ends up getting lumped with *N. faginata* (the ASV was named as *Microcera lavarum*, see table `LULU_which_taxa_reassigned.xlsx` for a breakdown of which taxa ASVs were reassigned to at different seq. sim. levels.) This ASV occurs in only one sample, but based on the sequence similarity is less likely an error.


### The need/desire for LULU was decided by first looking at sequence similarity clustering of ASVs within named taxa (i.e., the following routine was performed first)

Sequence similarity clustering of ASV seqs
```
#vsearch is loaded with itsxpress
conda activate itsxpress_3p
mkdir seq_sim_ASV_cluster
```
`vsearch --cluster_size` uses abundance based sorting
```
ID_vals=(0.95 0.96 0.97 0.98 0.99 1)
for i in ${ID_vals[@]}
do(
vsearch --cluster_size dada2_out/ASVs.fa \
--id $i \
--strand plus \
--sizein \
--sizeout \
--fasta_width 0 \
--uc seq_sim_ASV_cluster/ASVs.${i}.uc \
--centroids seq_sim_ASV_cluster/ASVs.centroids.${i}.fasta
)
done

conda deactivate
```
To filter to number of ASVs per taxon first extract list of ASVs that formed centroids, then use this list to filter taxa in R and summarize table as before
```
for i in ${ID_vals[@]}
do(
perl -ne 'print "$1\n" if/>(ASV_\d+);/' seq_sim_ASV_cluster/ASVs.centroids.$i.fasta > seq_sim_ASV_cluster/ASVs.centroids.$i.txt
echo $i
wc -l seq_sim_ASV_cluster/ASVs.centroids.$i.txt
)
done
```

```
Rscript ~/repo/neonectria_barcoding_012220/ecol/lowest_reliable_taxon_table_w_otus.r
```

