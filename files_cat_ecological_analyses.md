

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



