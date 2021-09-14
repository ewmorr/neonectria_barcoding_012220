
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")

Nf_v_Nd.bin.metadata.gt1K = filter(Nf_v_Nd.bin.metadata, total_seqs > 999 & sample != "BP79")


Nd.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)/nrow(Nf_v_Nd.bin.metadata.gt1K)
Nf.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)/nrow(Nf_v_Nd.bin.metadata.gt1K)
both.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "both") %>% summarize(n()) )/nrow(Nf_v_Nd.bin.metadata.gt1K)
both_when_one_present.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "both") %>% summarize(n()) )/( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence != "none") %>% summarize(n()))
none.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "none") %>% summarize(n()) )/nrow(Nf_v_Nd.bin.metadata.gt1K)

Nd.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)
Nf.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)
min_cooccur = Nf.num + Nd.num - nrow(Nf_v_Nd.bin.metadata.gt1K)
max_cooccur = min(Nd.num, Nf.num)

runif(171, 0, 1)
Nf_Nd_sxs = t(data.frame(Nf = Nf_v_Nd.bin.metadata.gt1K$Nf, Nd = Nf_v_Nd.bin.metadata.gt1K$Nd,
rand1 = round(runif(171, 0, 1)),
rand2 = round(runif(171, 0, 1)),
rand3 = round(runif(171, 0, 1))
))

colnames(Nf_Nd_sxs) = Nf_v_Nd.bin.metadata.gt1K$sample
mat_type = "spp_site"

Nf_Nd.cooccur.plug = cooccur(mat = Nf_Nd_sxs,
    type = mat_type,
    prob = "comb",
    thresh=TRUE,
    spp_names=T
)

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

Nf_v_Nd.bin.metadata.gt1K = filter(Nf_v_Nd.bin.metadata, total_seqs > 999)


Nd.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)/nrow(Nf_v_Nd.bin.metadata.gt1K)
Nf.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)/nrow(Nf_v_Nd.bin.metadata.gt1K)
both.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "both") %>% summarize(n()) )/nrow(Nf_v_Nd.bin.metadata.gt1K)
both_when_one_present.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "both") %>% summarize(n()) )/( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence != "none") %>% summarize(n()))
none.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "none") %>% summarize(n()) )/nrow(Nf_v_Nd.bin.metadata.gt1K)

Nd.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)
Nf.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)
min_cooccur = Nf.num + Nd.num - nrow(Nf_v_Nd.bin.metadata.gt1K)
max_cooccur = min(Nd.num, Nf.num)

runif(102, 0, 1)
Nf_Nd_sxs = t(data.frame(Nf = Nf_v_Nd.bin.metadata.gt1K$Nf, Nd = Nf_v_Nd.bin.metadata.gt1K$Nd,
rand1 = round(runif(102, 0, 1)),
rand2 = round(runif(102, 0, 1)),
rand3 = round(runif(102, 0, 1))
))

colnames(Nf_Nd_sxs) = Nf_v_Nd.bin.metadata.gt1K$sample
mat_type = "spp_site"

Nf_Nd.cooccur.tree = cooccur(mat = Nf_Nd_sxs,
type = mat_type,
prob = "comb",
thresh=TRUE,
spp_names=T
)


#No WI

Nf_v_Nd.bin.metadata.gt1K = filter(Nf_v_Nd.bin.metadata, total_seqs > 999 & Site != "WF1")


Nd.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)/nrow(Nf_v_Nd.bin.metadata.gt1K)
Nf.freq = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)/nrow(Nf_v_Nd.bin.metadata.gt1K)
both.freq = ( Nf_v_Nd.bin.metadata.gt1K %>% filter(occurence == "both") %>% summarize(n()) )/nrow(Nf_v_Nd.bin.metadata.gt1K)

Nd.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nd)
Nf.num = sum(Nf_v_Nd.bin.metadata.gt1K$Nf)
min_cooccur = Nf.num + Nd.num - nrow(Nf_v_Nd.bin.metadata.gt1K)
max_cooccur = min(Nd.num, Nf.num)

runif(81, 0, 1)
Nf_Nd_sxs = t(data.frame(Nf = Nf_v_Nd.bin.metadata.gt1K$Nf, Nd = Nf_v_Nd.bin.metadata.gt1K$Nd,
rand1 = round(runif(81, 0, 1)),
rand2 = round(runif(81, 0, 1)),
rand3 = round(runif(81, 0, 1))
))

colnames(Nf_Nd_sxs) = Nf_v_Nd.bin.metadata.gt1K$sample
mat_type = "spp_site"

Nf_Nd.cooccur.tree.noWF1 = cooccur(mat = Nf_Nd_sxs,
type = mat_type,
prob = "comb",
thresh=TRUE,
spp_names=T
)
