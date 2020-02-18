require(dada2)

seqtab.nochim = readRDS("intermediate_RDS/dada2_seq_table_no_chim.rds")
unite.ref <- "~/sh_general_release_s_all_04.02.2020/sh_general_release_dynamic_s_all_04.02.2020.fasta"
taxa.w_bootstraps <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = 24, tryRC = TRUE, outputBootstraps = T) #With latest R and tidyverse update this step can result in segfault and memalloc errors with tidyverse loaded. Need to load tidyverse after this step if needed (currently using for plotting and table manipulations

saveRDS(taxa.w_bootstraps, "intermediate_RDS/taxa_w_bootstraps.rds")
