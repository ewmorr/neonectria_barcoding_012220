require(dada2)

seqtab.nochim = readRDS("dada2_seq_table_no_chim.rds")
unite.ref <- "/Users/ericmorrison/blast_dbs/sh_general_release_dynamic_02.02.2019/sh_general_release_dynamic_02.02.2019.fasta"
taxa.w_bootstraps <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, outputBootstraps = T) #With latest R and tidyverse update this step can result in segfault and memalloc errors with tidyverse loaded. Need to load tidyverse after this step if needed (currently using for plotting and table manipulations

saveRDS(taxa.w_bootstraps, "intermediate_RDS/taxa_w_bootstraps.rds")
