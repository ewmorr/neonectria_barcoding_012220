

Data is housed in `~/GARNAS_neonectria_barcoding_files_cat_03242020`.
This is data that has been processed by first quality filtering and ITS extaction of individual .fastq files, then concatenating files from the same sample and running through dada2 algorithm (`sequence_processing_DADA2_README.md`). Finally, ASVs were filtered using LULU algorithm (`LULU_README.md`).
```
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
```

Make table of site-wide averages for ancillary data
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/site_averages_on_transects_data.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/sites_climate_dat.r
```
Analysis of samples and sequences dropped/retained at different rarefaction levels
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/sequence_subsampling_samples_seqs_dropped.r
```
Run NMDS at rarefaction levels of 1K and 5K, also save rarefied tables
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/run_and_save_NMDS.r
```
prelim figs (NMDS, site frequency Nf/Nd, richness etc)
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/prelim_figs_data_explore_files_cat.r
```
Table of lowest informative taxon
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/lowest_reliable_taxon_table.r
```


##### extract nectriaceae and/or neonectria ASV seqs and perform phylogenetic analyses
