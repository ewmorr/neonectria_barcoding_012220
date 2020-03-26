

Data is housed in `~/GARNAS_neonectria_barcoding_files_cat_03242020`
```
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
```

Make table of site-wide averages for ancillary data
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/site_averages_on_transects_data.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/sites_climate_dat.r
# The above script will calculate Tmax, MAT, Tmin, and ppt for sites
# then add to read_ASV_dat.r script (can then caluclate site-wide indicators (such as Nf frequency) as a function of these)
```
prelim figs (NMDS, site frequency Nf/Nd, richness etc)
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/prelim_figs_data_explore_files_cat.r
```
Table of lowest informative taxon
```
Rscript ~/repo/neonectria_barcoding_012220/files_cat/lowest_reliable_taxon_table.r
```



