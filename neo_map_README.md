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

cat neo_map/Nf.fa neo_map/Nd.fa > neo_map/Nf_Nd.fa
```
Merge quality filtered and ITS extracted sequences with vsearch and dereplicate. then map. Merging is performed with default minimum overlap and maxdifs.
```
mkdir neo_map/merged_seqs

#activate conda env with VSEARCH installed
conda activate itsxpress_3p

for i in run1_run2_files_cat_seqs/*R1.fastq.gz
do(
    sample_id=${i#*/}
    sample_id=${sample_id%_R1.fastq.gz}
    
    vsearch --fastq_mergepairs run1_run2_files_cat_seqs/${sample_id}_R1.fastq.gz \
        --threads 4 \
        --reverse run1_run2_files_cat_seqs/${sample_id}_R2.fastq.gz \
        --fastaout neo_map/merged_seqs/${sample_id}.fasta
)
done
```
Derpelicate sequences with "cluster" size output, and also label sequences with sample IDs
```
for i in neo_map/merged_seqs/*.fasta
do(
    sample_id=${i#*/}
    sample_id=${sample_id%.fasta}

    vsearch --derep_fulllength neo_map/merged_seqs/${sample_id}.fasta \
        --sizeout \
        --output neo_map/merged_seqs/${sample_id}.derep.fasta

    sed -i '' "s/^>.*/&;sample=${sample_id};/g" neo_map/merged_seqs/${sample_id}.derep.fasta
)
done
```
Cat merged, derepped, labeled sequences to single file for mapping
```
cat neo_map/merged_seqs/*.derep.fasta > neo_map/all_seqs.derep.fa
```
Use `vsearch --usearch_global` mapping to construct otutab
```
vsearch --usearch_global neo_map/all_seqs.derep.fa \
    --db neo_map/Nf_Nd.fa \
    --sizein \
    --id 0.97 \
    --alnout neo_map/all_seqs.derep.aln.txt \
    --otutabout neo_map/all_seqs.derep.otu_tab.txt

sed -i '' 's/#OTU ID/OTUID/' neo_map/all_seqs.derep.otu_tab.txt

conda deactivate
```
