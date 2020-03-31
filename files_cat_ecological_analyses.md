

Data is housed in `~/GARNAS_neonectria_barcoding_files_cat_03242020`
```
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
```

Make table of site-wide averages for ancillary data
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/site_averages_on_transects_data.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/sites_climate_dat.r
```
Run NMDS at rarefaction levels of 1K and 5K, also save rarefied tables
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/run_and_save_NMDS.r
```
prelim figs (NMDS, site frequency Nf/Nd, richness etc)
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/prelim_figs_data_explore_files_cat.r
```
Table of lowest informative taxon
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/lowest_reliable_taxon_table.r
```
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
Rscript ~/repo/neonectria_barcoding_012220/files_cat/lowest_reliable_taxon_table_w_otus.r
```

##### extract nectriaceae and/or neonectria ASV seqs and perform phylogenetic analyses 
