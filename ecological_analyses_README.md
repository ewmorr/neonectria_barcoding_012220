

Data is housed in `~/GARNAS_neonectria_barcoding_files_cat_03242020`.
This is data that has been processed by first quality filtering and ITS extaction of individual .fastq files, then concatenating files from the same sample and running through dada2 algorithm (`sequence_processing_DADA2_README.md`). Finally, ASVs were filtered using LULU algorithm (`LULU_README.md`).
```
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
```
#### Site level data
Make table of site-wide averages for ancillary data
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/site_averages_on_transects_data.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/sites_climate_dat.r
```
Calculate GDD based seasonal metrics using PRISM daily data. Some of the commands in the script for downloading dailys are hard-c-ded and should be checked before running.
```
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/tmin_tmax_dailys_site_extract.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/tmin_tmax_dailys_GDD_calcs.r
Rscript ~/repo/neonectria_barcoding_012220/PRISM_analysis/long-term_and_short-term_GDD_avgs.r
```
Pairwise comps of site dat
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/pairwise_comps_site_dat.r
```
#### Analysis of samples and sequences dropped/retained at different rarefaction levels
```
mkdir rarefaction_figs
Rscript ~/repo/neonectria_barcoding_012220/ecol/sequence_subsampling_samples_seqs_dropped.r
```
#### Run NMDS at rarefaction levels of 1K and 5K, also save rarefied tables. This is performed after filtering out ASVs that occur in only one sample (`read_ASV_dat.LULU_tab.r`)
First running stress v k
```
mkdir NMDS_fits
Rscript ~/repo/neonectria_barcoding_012220/ecol/NMDS_stress_v_k.r
```
No convergence at k=2, filtering ASVs at higher min frequency (e.g., 2 sample min occurence, 3 sample min occurence) results in no reduction in stress values at a given k. using k = 3 and ecluding singletons (single sample) ASVs
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/run_and_save_NMDS.r
```
#### prelim figs
(NMDS, site frequency Nf/Nd, richness etc)
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/prelim_figs_data_explore_files_cat.r
```
Table of lowest informative taxon
```
Rscript ~/repo/neonectria_barcoding_012220/ecol/lowest_reliable_taxon_table.r
```
GAM fits and plots
```
mkdir GAM_fits
Rscript ~/repo/neonectria_barcoding_012220/ecol/GAM_fits_NMDS.r
```

### sum samples at tree level and rerunning above analyses (rarefaction, NMDS, GAM fits)
Rarefaction
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/sequence_subsampling_samples_seqs_dropped.r
```
NMDS stress v k
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/NMDS_stress_v_k.r
```
Run NMDS, k = 3, singletons removed
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/run_and_save_NMDS.r
```
Prelim figs
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/prelim_figs_data_explore_sum_trees.r
```
GAM fits
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/GAM_fits_NMDS.site_means_cov.tree_level_comm.r
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/GAM_fits_NMDS.tree_level_cov.tree_level_comm.r
```
Added GAM fits for PRISM daily GDD calcs
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/GAM_fits_NMDS.site_GDDclim_cov.tree_level_comm.r
```
#### effects of number of plugs on richness and comm. comp.

```
mkdir n_plug_effects
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/n_plug_effects.r
```

#### Site level comm comp for landscape/continental vars
GAM fits on site level dat with community aggregated by site
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/run_and_save_NMDS.site_sum.r
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/GAM_fits_NMDS.site_means_cov.site_agg_comm.r
```
Mantel correlogram
```
Rscript ~/repo/neonectria_barcoding_012220/sum_trees/mantel_correlog.site_sum.r
```



##### extract nectriaceae and/or neonectria ASV seqs and perform phylogenetic analyses
