
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")

Nf_v_Nd.bin.metadata.gt1K = filter(Nf_v_Nd.bin.metadata, total_seqs > 999)

#Remove NA row for NeoCoverage
Nf_v_Nd.bin.metadata.gt1K = Nf_v_Nd.bin.metadata.gt1K[-87,]

Nf_v_Nd.bin.metadata.gt1K$NeoCoverage[Nf_v_Nd.bin.metadata.gt1K$NeoCoverage > 0] = 1

plot(Nf ~ NeoCoverage, data = Nf_v_Nd.bin.metadata.gt1K)

plot(Nd ~ NeoCoverage, data = Nf_v_Nd.bin.metadata.gt1K)

Nf_Nd_peri = data.frame(
    Nf = Nf_v_Nd.bin.metadata.gt1K$Nf,
    Nd = Nf_v_Nd.bin.metadata.gt1K$Nd,
    occurrence = Nf_v_Nd.bin.metadata.gt1K$occurence,
    peri = Nf_v_Nd.bin.metadata.gt1K$NeoCoverage
)

Nf_Nd_peri %>% group_by(Nf, peri) %>% summarize(n = n())

Nf_Nd_peri %>% group_by(Nd, peri) %>% summarize(n = n())

occ_by_peri = Nf_Nd_peri %>% group_by(occurrence, peri) %>% summarize(n = n())

Nf_Nd_peri %>% group_by(peri) %>% summarize(n = n())

write.table(occ_by_peri, "occurence_w_plug_peri.txt", sep = "\t", quote = F, row.names = F)
