Mapping of original reads after quality filtering and ITS extraction against *Neonectria* representative sequences. The representative sequences are the result of DADA2 before LULU filtering in order to provide the most complete set of possible sequences, and identified as *Neonectria* spp by  RDP naive bayesian taxonomic assignment against UNITE db. Mapping performed using VSEARCH  `--usearch_global`
```
cd ~/GARNAS_neonectria_barcoding_files_cat_03242020
mkdir neo_map
```
Get IDs of N. faginata and N. ditissima, and extract from the ASV sequence file
```
grep "s__faginata" dada2_out/ASVs_taxonomy.tsv | cut -f 1 > neo_map/Nf_ids.txt
grep "s__ditissima" dada2_out/ASVs_taxonomy.tsv | cut -f 1 > neo_map/Nd_ids.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' neo_map/Nf_ids.txt dada2_out/ASVs.fa > neo_map/Nf.fa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' neo_map/Nd_ids.txt dada2_out/ASVs.fa > neo_map/Nd.fa
```
Merge quality filtered and ITS extracted sequences with vsearch, then map. merging is performed with default minimum overlap and maxdifs
```
mkdir merged_seqs
mkdir neo_map/vsearch_aln_Nd
mkdir neo_map/vsearch_aln_Nf

#actiate conda env with VSEARCH installed
conda activate itsxpress_3p

for i in run1_run2_files_cat_seqs/*R1.fastq.gz
do(
    sample_id=${i#*/}
    sample_id=${sample_id%_R1.fastq.gz}
    vsearch --fastq_mergepairs run1_run2_files_cat_seqs/${sample_id}_R1.fastq.gz \
        --threads 4 \
        --reverse run1_run2_files_cat_seqs/${sample_id}_R2.fastq.gz \
        --fastaout merged_seqs/${sample_id}.fasta
)
done

#Run Nf and Nd separate
for i in merged_seqs/*
do(
    sample_id=${i##*/}
    sample_id=${sample_id%.fasta}
    vsearch --usearch_global $i --db neo_map/Nf.fa --id 0.97 --alnout neo_map/vsearch_aln_Nf/$sample_id.Nf.aln.txt
    vsearch --usearch_global $i --db neo_map/Nd.fa --id 0.97 --alnout neo_map/vsearch_aln_Nd/$sample_id.Nd.aln.txt

    echo $sample_id
)
done
```
Make OTU tables based on number of hits to database seqs
```
perl ~/repo/neonectria_barcoding_012220/neo_map/count_mapping_hits.pl neo_map/vsearch_aln_Nd
perl ~/repo/neonectria_barcoding_012220/neo_map/count_mapping_hits.pl neo_map/vsearch_aln_Nf

