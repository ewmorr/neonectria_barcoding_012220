library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)

source("~/ggplot_theme.txt")

fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    #reformat zeros
    l <- gsub("0e\\+00","0",l)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
}

get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

setwd("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320")
system("mkdir separate_denoising_per_run")

seqtab.nochim.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/taxa_w_bootstraps.rds")

seqtab.nochim.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/taxa_w_bootstraps.rds")

#set sample names to match metadata
rownames(seqtab.nochim.r1) = unname(sapply(rownames(seqtab.nochim.r1), get.sample.name))
rownames(seqtab.nochim.r1) = paste0("runOne", rownames(seqtab.nochim.r1))
rownames(seqtab.nochim.r2) = unname(sapply(rownames(seqtab.nochim.r2), get.sample.name))
rownames(seqtab.nochim.r2) = paste0("runTwo", rownames(seqtab.nochim.r2))

#Read metadata
id_bench_map = read.table("sample_data/sample_mapping.txt", header = T)
metadata_map = read.table("sample_data/metadata.txt", header = T)
survey_dat = read.table("sample_data/trees_site_survey_data.txt", header = T, sep = "\t")
neo_cov = read.table("sample_data/plug_neonectria_coverage.txt", header = T)
site_info = read.csv("sample_data/site_info.csv")

colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>% left_join(., site_info, by = "Site")

full_metadata = full_join(metadata_ordered, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))

#############
#join tables#

#full count tab
seqtab.nochim.join_by_ASVseq = full_join(
    data.frame(asv.seq = rownames(t(seqtab.nochim.r1)), t(seqtab.nochim.r1)),
    data.frame(asv.seq = rownames(t(seqtab.nochim.r2)), t(seqtab.nochim.r2)),
    by = "asv.seq"
)
rownames(seqtab.nochim.join_by_ASVseq) = seqtab.nochim.join_by_ASVseq$asv.seq
seqtab.nochim.join_by_ASVseq$asv.seq = NULL
seqtab.nochim.join_by_ASVseq[is.na(seqtab.nochim.join_by_ASVseq)] = 0

#full taxa table
#Its seems after doing some joins that some taxonomic names are not matching between different seqs, and this may be attribued at the Kingdom level (based on error from full_join). Will need to explore this

asv_tax.join_by_ASVseq = rbind(
    data.frame(asv.seq = rownames(taxa.w_bootstraps.r1$tax), taxa.w_bootstraps.r1$tax),
    data.frame(asv.seq = rownames(taxa.w_bootstraps.r2$tax), taxa.w_bootstraps.r2$tax)
) %>% unique

seqtab.nochim.join_by_ASVseq.w_tax = full_join(
    data.frame(asv.seq = rownames(seqtab.nochim.join_by_ASVseq), seqtab.nochim.join_by_ASVseq),
    data.frame(asv.seq = rownames(asv_tax.join_by_ASVseq), asv_tax.join_by_ASVseq),
    by = "asv.seq"
)


#asv_tax.join_by_ASVseq = full_join(
#data.frame(asv.seq = rownames(taxa.w_bootstraps.r1$tax), taxa.w_bootstraps.r1$tax),
#data.frame(asv.seq = rownames(taxa.w_bootstraps.r2$tax), taxa.w_bootstraps.r2$tax)
#, by = "asv.seq"
#)

#full_join(
#data.frame(asv.seq = rownames(taxa.w_bootstraps.r1$tax), taxa.w_bootstraps.r1$tax),
#data.frame(asv.seq = rownames(taxa.w_bootstraps.r2$tax), taxa.w_bootstraps.r2$tax)
#) %>% as_tibble
#asv_tax.join_by_ASVseq %>% filter(Species.x != Species.y)

asv_tax.join_by_ASVseq.Nf = filter(asv_tax.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__faginata")
asv_tax.join_by_ASVseq.Nf$asv.seq %>% length
asv_tax.join_by_ASVseq.Nf$asv.seq %>% unique %>% length

asv_tax.join_by_ASVseq.Nd = filter(asv_tax.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__ditissima")
asv_tax.join_by_ASVseq.Nd$asv.seq %>% length
asv_tax.join_by_ASVseq.Nd$asv.seq %>% unique %>% length

######################################################################
#tables of run one vs. run two values by sample (i.e. metadata.label)#

#total sequence counts after full dada2
sample_counts = data.frame(sample = colnames(seqtab.nochim.join_by_ASVseq), count = colSums(seqtab.nochim.join_by_ASVseq))
sample_counts.meta = left_join(
    full_metadata %>%
        filter(seq.rep == "n" & locus == "ITS2") %>%
        select(c("metadata.label", "sample", "run.seq")) ,
    sample_counts,
    by = "sample"
)
sample_counts.meta[is.na(sample_counts.meta)] = 0

sample_counts.meta.wide_run = pivot_wider(sample_counts.meta, id_cols = metadata.label, names_from = run.seq, values_from = count)

qqnorm(residuals(lm(log10(two+1) ~ log10(one+1), data = sample_counts.meta.wide_run)))
plot(residuals(lm(log10(two+1) ~ log10(one+1), data = sample_counts.meta.wide_run)))
summary(lm(log10(two+1) ~ log10(one+1), data = sample_counts.meta.wide_run))

fit = lm(log10(two+1) ~ log10(one+1), data = sample_counts.meta.wide_run)

p = ggplot(sample_counts.meta.wide_run, aes(x = one + 1, y = two + 1)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific, limits = c(-1,110000)) +
scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "run one sequence count", y = "run two sequence count",
title = paste(
"Sequence count by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
) +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_seq_count.pdf", width = 8, height = 6)
print(p)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(seqtab.nochim.join_by_ASVseq), seqtab.nochim.join_by_ASVseq) %>%
    filter(asv.seq %in% asv_tax.join_by_ASVseq.Nf$asv.seq)
seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(seqtab.nochim.join_by_ASVseq),
richness = colSums(seqtab.nochim.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.join_by_ASVseq)
)

Nf_count.meta = left_join(
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
select(c("metadata.label", "sample", "run.seq")) ,
Nf_count,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(count, richness, Nf))

p = ggplot(Nf_count.meta.wide %>% filter(count_one != 0 & count_two != 0), aes(Nf_one/count_one, Nf_two/count_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_log10(labels = fancy_scientific) +
#scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "run one N. faginata rel. abd.", y = "run two N. faginata rel. abd.") +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_Nf_RA.pdf", width = 8, height = 6)
print(p)
dev.off()


#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(seqtab.nochim.join_by_ASVseq), seqtab.nochim.join_by_ASVseq) %>%
filter(asv.seq %in% asv_tax.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(seqtab.nochim.join_by_ASVseq),
richness = colSums(seqtab.nochim.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.join_by_ASVseq)
)

Nd_count.meta = left_join(
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq")) ,
Nd_count,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(count, richness, Nd))

p = ggplot(Nd_count.meta.wide %>% filter(count_one != 0 & count_two != 0), aes(Nd_one/count_one, Nd_two/count_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
scale_x_sqrt(labels = fancy_scientific) +
scale_y_sqrt(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "run one N. ditissima rel. abd.", y = "run two N. ditissima rel. abd.") +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_Nd_RA.pdf", width = 8, height = 6)
print(p)
dev.off()




########################################################
#Pairwise asv_tabs for rarefying samples by run minimum#
########################################################
require(vegan)

#table for selecting sample pairs
sample_map.wide_run = full_metadata %>% filter(!is.na(metadata.label) & seq.rep == "n" & locus == "ITS2") %>% select(c(sample, run.seq, metadata.label))  %>%
    pivot_wider(., id_cols = metadata.label, names_from = run.seq, values_from = sample)

sample_map.wide_run = sample_map.wide_run %>% filter(!is.na(one) & !is.na(two) )

seqtab.nochim.join_by_ASVseq.t = t(seqtab.nochim.join_by_ASVseq)
seqtab.nochim.join_by_ASVseq.t.sample = data.frame(sample = rownames(seqtab.nochim.join_by_ASVseq.t), seqtab.nochim.join_by_ASVseq.t)

#empty matrix with rows as asvs and cols as samples
seqtab.nochim.sample_rare_by_run_min = matrix(
    nrow = length(rownames(seqtab.nochim.join_by_ASVseq.t)),
    ncol = length(colnames(seqtab.nochim.join_by_ASVseq.t)),
    dimnames = list(rownames(seqtab.nochim.join_by_ASVseq.t), colnames(seqtab.nochim.join_by_ASVseq.t))
)

#simplest method is to call row by var name (bc dplyr filter NSE variable interpolation)
#This will break if the rownames in the empty matrix are not properly names, so need to run the above call if the matrix is modified
for(i in 1:length(sample_map.wide_run$metadata.label)){
    #skip sample pairs where one is missing
    if(any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$one[i]) == FALSE |
    any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$two[i]) == FALSE){
        next
    }
    
    asv_tab.temp.runOne = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$one[i]), ]
    asv_tab.temp.runTwo = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$two[i]), ]
    min_seqs = min(sum(asv_tab.temp.runOne), sum(asv_tab.temp.runTwo))
    
    asv_tab.temp = rbind(asv_tab.temp.runOne, asv_tab.temp.runTwo)
    rownames(asv_tab.temp) = c(as.character(sample_map.wide_run$one[i]), as.character(sample_map.wide_run$two[i]))
    asv_tab.temp.rare = rrarefy(asv_tab.temp, min_seqs)
    seqtab.nochim.sample_rare_by_run_min[c(as.character(sample_map.wide_run$one[i]), as.character(sample_map.wide_run$two[i])),] =
        asv_tab.temp.rare[c(as.character(sample_map.wide_run$one[i]), as.character(sample_map.wide_run$two[i])),]
}

seqtab.nochim.sample_rare_by_run_min.rmNA = drop_na(as.data.frame(seqtab.nochim.sample_rare_by_run_min))


#####################################################
#Calculate vals and create pairwise DFs for rarefied#

######################
#rarefied Simoson's D#
seqtab.nochim.sample_rare_by_run_min.rmNA.simp = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "simpson")


asv_simpson = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
simp.D = unname(1 - seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_simpson.meta = left_join(asv_simpson,
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq")),
by = "sample"
)

asv_simpson.meta.wide = pivot_wider(asv_simpson.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(simp.D, count))

qqnorm(residuals(lm(simp.D_two ~ simp.D_one, data = asv_simpson.meta.wide %>% filter(count_one >= 100))))
plot(residuals(lm(simp.D_two ~ simp.D_one, data = asv_simpson.meta.wide %>% filter(count_one >= 100))))
summary(lm(simp.D_two ~ simp.D_one, data = asv_simpson.meta.wide %>% filter(count_one >= 100)))

fit = lm(simp.D_two ~ simp.D_one, data = asv_simpson.meta.wide %>% filter(count_one >= 100) )

p = ggplot(asv_simpson.meta.wide  %>% filter(count_one >= 100), aes(simp.D_one, simp.D_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "run one Simpson D", y = "run two Simpson D",
    title = paste(
        "Simpson D by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
        " Slope = ",signif(fit$coef[[2]], 2),
        " P =",signif(summary(fit)$coef[2,4], 3)
    )
)

pdf("separate_denoising_per_run/run1_v_run2_Simpson_D.pdf", width = 8, height = 6)
print(p)
dev.off()

######################
#rarefied Shannon H#
seqtab.nochim.sample_rare_by_run_min.rmNA.shan = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "shannon")


asv_shannon = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
shannon = unname(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_shannon.meta = left_join(asv_shannon,
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq")),
by = "sample"
)

asv_shannon.meta.wide = pivot_wider(asv_shannon.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(shannon, count))

qqnorm(residuals(lm(shannon_two ~ shannon_one, data = asv_shannon.meta.wide %>% filter(count_one >= 100))))
plot(residuals(lm(shannon_two ~ shannon_one, data = asv_shannon.meta.wide %>% filter(count_one >= 100))))
summary(lm(shannon_two ~ shannon_one, data = asv_shannon.meta.wide %>% filter(count_one >= 100)))

fit = lm(shannon_two ~ shannon_one, data = asv_shannon.meta.wide %>% filter(count_one >= 100) )

p = ggplot(asv_shannon.meta.wide %>% filter(count_one >= 100), aes(shannon_one, shannon_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "run one Shannon H", y = "run two Shannon H",
title = paste(
    "Shannon H by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
    " Slope = ",signif(fit$coef[[2]], 2),
    " P =",signif(summary(fit)$coef[2,4], 3)
    )
)


pdf("separate_denoising_per_run/run1_v_run2_shannon.pdf", width = 8, height = 6)
print(p)
dev.off()


##############
#asv richness#
asv_richness = data.frame(
sample = rownames(seqtab.nochim.sample_rare_by_run_min.rmNA),
richness = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA > 0),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_richness.meta = left_join(asv_richness,
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
select(c("metadata.label", "sample", "run.seq")) ,
by = "sample"
)
asv_richness.meta[is.na(asv_richness.meta)] = 0

asv_richness.meta.wide = pivot_wider(asv_richness.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(count, richness))

qqnorm(residuals(lm(richness_two ~ richness_one, data = asv_richness.meta.wide %>% filter(count_one >= 100 & count_two >= 100))))
plot(residuals(lm(richness_two ~ richness_one, data = asv_richness.meta.wide %>% filter(count_one >= 100 & count_two >= 100))))
summary(lm(richness_two ~ richness_one, data = asv_richness.meta.wide %>% filter(count_one >= 100 & count_two >= 100)))

fit = lm(richness_two ~ richness_one, data = asv_richness.meta.wide %>% filter(count_one >= 100 & count_two >= 100))

require(gridExtra)

p1 = ggplot(asv_richness.meta.wide %>% filter(count_one >= 100 & count_two >= 100), aes(richness_one, richness_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "run one ASV richness", y = "run two ASV richness",
    title = paste(
        "ASV richness by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
        " Slope = ",signif(fit$coef[[2]], 2),
        " P =",signif(summary(fit)$coef[2,4], 3)
    )
)

p2 = ggplot(asv_richness.meta.wide %>% filter(count_one != 0 & count_two != 0), aes(richness_one, richness_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "run one ASV richness log10+1", y = "run two ASV richness log10+1", title = "ASV richness by sample log10") +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_richness.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.join_by_ASVseq.Nf$asv.seq)
seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nf_count.meta = left_join(
Nf_count,
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq")) ,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(count, richness, Nf))

qqnorm(residuals(lm(log10(Nf_two+1) ~ log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100))))
plot(residuals(lm(log10(Nf_two+1) ~ log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100))))
summary(lm(log10(Nf_two+1) ~ log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100)))

qqnorm(residuals(lm(Nf_two/count_two ~ I(Nf_one/count_one), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100)))
summary(lm(Nf_two/count_two ~ I(Nf_one/count_one), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100)))

#proportional counts
Nf.binomial_glm = glm(Nf_two/count_two ~ I(Nf_one/count_one*100), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), family = binomial)
1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance
)
1 - pchisq(summary(Nf.binomial_glm)$deviance,
summary(Nf.binomial_glm)$df.residual
)

#predict(Nf.binomial_glm, )
#models for raw counts (but proportional is probably fine/most appropriate
#Nf.poisson_glm = glm(Nf_two ~ log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), family = poisson)

#plot(residuals(Nf.poisson_glm))
#qqnorm(residuals(Nf.poisson_glm))


#require(MASS)
#Nf.poisson_glm.nb = glm.nb(Nf_two ~ log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100))
#summary(Nf.poisson_glm.nb)

#1 - pchisq(summary(Nf.poisson_glm.nb)$deviance,
#summary(Nf.poisson_glm.nb)$df.residual
#)

#Nf.zeroinfl = zeroinfl(Nf_two ~ log10(Nf_one+1)|log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), dist = "poisson", link = "logit")
#Nf.zeroinfl.nb = zeroinfl(Nf_two ~ log10(Nf_one+1)|log10(Nf_one+1), data = Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), dist = "negbin")
#summary(Nf.zeroinfl.nb)

p1 = ggplot(Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), aes(Nf_one/count_one, Nf_two/count_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
my_gg_theme +
labs(x = "run one N. faginata rel. abd.", y = "run two N. faginata rel. abd.",
title = paste(
"N. faginata rel abd by sample\nR2 = ", signif(1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nf.binomial_glm)$coef[2,4], 3)
)
)

p2 = ggplot(Nf_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), aes(Nf_one+1, Nf_two+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "run one N. faginata count log10+1", y = "run two N. faginata count log10+1", title = "N. faginata counts") +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_Nf_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nd_count.meta = left_join(
Nd_count,
full_metadata %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq")) ,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = metadata.label, names_from = run.seq, values_from = c(count, richness, Nd))

qqnorm(residuals(lm(log10(Nd_two+1) ~ log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100 & Nd_one+Nd_two != 0))))
plot(residuals(lm(log10(Nd_two+1) ~ log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100 & Nd_one+Nd_two != 0))))
summary(lm(log10(Nd_two+1) ~ log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100 & Nd_one+Nd_two != 0)))
summary(lm(Nd_two/count_two ~ I(Nd_one/count_one), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100)))

#proportional counts
Nd.binomial_glm = glm(Nd_two/count_two ~ I(Nd_one/count_one), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), family = binomial)
1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance
)
1 - pchisq(summary(Nd.binomial_glm)$deviance,
summary(Nd.binomial_glm)$df.residual
)

#log 10 transform predictor counts for log-linearity with dependent
Nd.poisson_glm = glm(Nd_two ~ log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), family = poisson)

# if > 0.05 model is appropriate
1 - pchisq(summary(Nd.poisson_glm)$deviance,
summary(Nd.poisson_glm)$df.residual
)

Nd.poisson_glm.nb = glm.nb(Nd_two ~ log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100))
summary(Nd.poisson_glm.nb)

1 - pchisq(summary(Nd.poisson_glm.nb)$deviance,
summary(Nd.poisson_glm.nb)$df.residual
)

Nd.zeroinfl = zeroinfl(Nd_two ~ log10(Nd_one+1)|log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), dist = "poisson", link = "logit")

#The second predictor term allows different levels of zero inflation for different levels of predictor
Nd.zeroinfl.nb = zeroinfl(Nd_two ~ log10(Nd_one+1)|log10(Nd_one+1), data = Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), dist = "negbin")
summary(Nd.zeroinfl.nb)
#if theta is sig the nb model is appropriate, and vice versa

p1 = ggplot(Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), aes(Nd_one/count_one, Nd_two/count_two)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_sqrt(labels = fancy_scientific) +
scale_y_sqrt(labels = fancy_scientific) +
my_gg_theme +
labs(x = "run one N. ditissima rel. abd.", y = "run two N. ditissima rel. abd.",
title = paste(
"N. ditissima rel abd by sample\nR2 = ", signif(1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nd.binomial_glm)$coef[2,4], 1)
)
)

p2 = ggplot(Nd_count.meta.wide %>% filter(count_one >= 100 & count_two >= 100), aes(Nd_one+1, Nd_two+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "run one N. ditissima count log10+1", y = "run two N. ditissima count log10+1", title = "N. ditissima counts") +
my_gg_theme

pdf("separate_denoising_per_run/run1_v_run2_Nd_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()





##################################
#Prep OTU table factors for adonis
sample_factors = left_join(
data.frame(sample = rownames(seqtab.nochim.sample_rare_by_run_min.rmNA),
full_metadata,
by = "sample"
)



adonis(seqtab.nochim.sample_rare_by_run_min.rmNA ~ sample_factors$run.seq + sample_factors$metadata.label, strata = sample_factors$metadata.label)

#Call:
#adonis(formula = seqtab.nochim.sample_rare_by_run_min.rmNA ~      sample_factors$run.seq + sample_factors$metadata.label, strata = sample_factors$metadata.label)
#
#Blocks:  strata
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#sample_factors$run.seq          1     0.057 0.05663  1.5326 0.00036  0.002 **
#sample_factors$metadata.label 176   149.398 0.84885 22.9726 0.95794  0.002 **
#Residuals                     176     6.503 0.03695         0.04170
#Total                         353   155.958                 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(log10(seqtab.nochim.sample_rare_by_run_min.rmNA+1) ~ sample_factors$run.seq + sample_factors$metadata.label, strata = sample_factors$metadata.label)

#Call:
#adonis(formula = log10(seqtab.nochim.sample_rare_by_run_min.rmNA +      1) ~ sample_factors$run.seq + sample_factors$metadata.label,      strata = sample_factors$metadata.label)
#
#Blocks:  strata
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#sample_factors$run.seq          1     0.644 0.64363  10.358 0.00453  0.001 ***
#sample_factors$metadata.label 176   130.537 0.74169  11.935 0.91851  0.001 ***
#Residuals                     176    10.937 0.06214         0.07696
#Total                         353   142.118                 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(vegdist(seqtab.nochim.sample_rare_by_run_min.rmNA, method = "bray", binary = T) ~ sample_factors$run.seq + sample_factors$metadata.label, strata = sample_factors$metadata.label)

#Call:
#adonis(formula = vegdist(seqtab.nochim.sample_rare_by_run_min.rmNA,      method = "bray", binary = T) ~ sample_factors$run.seq + sample_factors$metadata.label,      strata = sample_factors$metadata.label)
#
#Blocks:  strata
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#sample_factors$run.seq          1     1.504 1.50393 14.6846 0.01071  0.001 ***
#sample_factors$metadata.label 176   120.852 0.68666  6.7046 0.86089  0.001 ***
#Residuals                     176    18.025 0.10242         0.12840
#Total                         353   140.381                 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(seqtab.nochim.sample_rare_by_run_min.rmNA ~ sample_factors$run.seq + sample_factors$metadata.label)
adonis(seqtab.nochim.sample_rare_by_run_min.rmNA ~ sample_factors$run.seq, strata = sample_factors$metadata.label)


###############################################
#Calculate Nf and Nd detection agreement between runs#

agree_or_not = function(x){
    if(x[1] > 0 & x[2] > 0){
        return(1)
    }else if(x[1] == 0 & x[2] == 0){
        return(1)
    }else{
        return(0)
    }
}

Nf_count.meta.wide$detect_agree = apply(data.frame(Nf_count.meta.wide$Nf_one, Nf_count.meta.wide$Nf_two), 1, agree_or_not)

summary(glm(detect_agree ~ log10(count_one), data = Nf_count.meta.wide, family = binomial))
summary(glm(detect_agree ~ log10(Nf_one+Nf_two+1), data = Nf_count.meta.wide, family = binomial))
summary(glm(detect_agree ~ I((Nf_one+Nf_two)/count_one), data = Nf_count.meta.wide, family = binomial))


p1 = ggplot(Nf_count.meta.wide %>% filter(Nf_one+Nf_two != 0), aes(count_one, detect_agree, color = (Nf_one+Nf_two))) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "sequence count", y = "detected in both runs", color = "run1 + run2\nNf count") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

p2 = ggplot(Nf_count.meta.wide %>% filter(Nf_one+Nf_two != 0), aes((Nf_one+Nf_two), detect_agree, color = count_one)) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "run1 + run2 Nf count", y = "detected in both runs", color = "sequence\ncount") +
my_gg_theme +
theme(legend.title = element_text(size = 18))


p4 = ggplot(Nf_count.meta.wide %>% filter(Nf_one+Nf_two != 0), aes((Nf_one+Nf_two)/(count_one+count_two), detect_agree, color = count_one)) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
#scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "avg Nf rel. abd.", y = "detected in both runs", color = "sequence\ncount") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

p3 = ggplot(Nf_count.meta.wide %>% filter(Nf_one+Nf_two != 0), aes(count_one, detect_agree, color = (Nf_one+Nf_two)/(count_one+count_two))) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
labs(x = "sequence count", y = "detected in both runs", color = "avg Nf\nrel. abd.") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

pdf("separate_denoising_per_run/Nf_detection_by_seq_counts.pdf", width = 16, height = 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

Nd_count.meta.wide$detect_agree = apply(data.frame(Nd_count.meta.wide$Nd_one, Nd_count.meta.wide$Nd_two), 1, agree_or_not)

summary(glm(detect_agree ~ log10(count_one), data = Nd_count.meta.wide, family = binomial))
summary(glm(detect_agree ~ log10(Nd_one+Nd_two+1), data = Nd_count.meta.wide, family = binomial))
summary(glm(detect_agree ~ I((Nd_one+Nd_two+1)/count_one), data = Nd_count.meta.wide, family = binomial))

p1 = ggplot(Nd_count.meta.wide %>% filter(Nd_one+Nd_two != 0), aes(count_one, detect_agree, color = (Nd_one+Nd_two))) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "sequence count", y = "detected in both runs", color = "run1 + run2\nNd count") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

p2 = ggplot(Nd_count.meta.wide %>% filter(Nd_one+Nd_two != 0), aes(Nd_one+Nd_two, detect_agree, color = count_one)) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "run1 + run2 Nd count", y = "detected in both runs", color = "sequence\ncount") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

p4 = ggplot(Nd_count.meta.wide %>% filter(Nd_one+Nd_two != 0), aes((Nd_one+Nd_two)/(count_one+count_two), detect_agree, color = count_one)) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
#scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "log10") +
labs(x = "avg Nd rel. abd.", y = "detected in both runs", color = "sequence\ncount") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

p3 = ggplot(Nd_count.meta.wide %>% filter(Nd_one+Nd_two != 0), aes(count_one, detect_agree, color = (Nd_one+Nd_two)/(count_one+count_two))) +
geom_point(size = 3, alpha = 0.5, position = position_jitter(height = 0.025)) +
#stat_smooth(method = "glm", method.args = list(family = "binomial")) +
scale_y_continuous(breaks = c(0,1)) +
scale_x_log10(labels = fancy_scientific, limits = c(0.75, 120000), breaks = c(1,10,100,1000,10^4,10^5)) +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
labs(x = "sequence count", y = "detected in both runs", color = "avg Nd\nrel. abd.") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

pdf("separate_denoising_per_run/Nd_detection_by_seq_counts.pdf", width = 16, height = 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


