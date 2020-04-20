
mkdirs and rename R1 and R2 (this performs a "reverse complement" in essence, so that ITSxpress works properly)
```
cd ~/GARNAS_neonectria_barcoding_012220


mkdir original_reads

for i in cobb.sr.unh.edu/managed/*P3/reads/Sample*/*fastq.gz;do( mv $i ./original_reads/;)done

mkdir R1_R2_switched

for i in original_reads/*R1*; do(DIR=${i%/*}; FILE=${i##*/}; BEFFILE=${FILE%R1*}; AFTFILE=${FILE##*R1}; cp $DIR/${BEFFILE}R1${AFTFILE} ./R1_R2_switched/${BEFFILE}R2${AFTFILE}; cp $DIR/${BEFFILE}R2${AFTFILE} ./R1_R2_switched/${BEFFILE}R1${AFTFILE}); done

```

```
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS
```
Run the pre-processing script.
```
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/rm_primers_and_qual_filter.r
```
Or open an R terminal and run the commands in `rm_primers_and_qual_filter.r` (This stuff should mostly run fine on a reasonably good laptop)

Next run itsxpress on the filtered reads. This should use the version of itsxpress that has been hacked to trim 3' read ends, other wise we end up with staggered read merges that have conserved 5.8S or LSU regions remaining (not good)

The `hmmsearch` step can take quite a while, so this is one place to use a cluster.

```
conda activate itsxpress_3p

cd ~/GARNAS_neonectria_barcoding_012220/
INDIR=R1_R2_switched/filtN/cutadapt/filtered
OUTDIR=R1_R2_switched/itsxpress
mkdir $OUTDIR

for i in $INDIR/*R1*
do(
    FILE=${i##*/}
    BEFFILE=${FILE%R1*}
    AFTFILE=${FILE##*R1}
    R1=$FILE
    R2=${BEFFILE}R2${AFTFILE}
    echo $R1
    if [ -f $OUTDIR/$R2 ]
    then
        continue
    fi

    itsxpress \
    --fastq $INDIR/$R1 --fastq2 $INDIR/$R2 \
    --outfile $OUTDIR/$R1 --outfile2 $OUTDIR/$R2 \
    --region ITS2 --taxa 'Fungi' --cluster_id 1 \
    --threads 4 \
    --log itsxpress.log
)
done
```
The stdout is sent to file bc I have the program set to output the trimmed read len of each fwd and reverse read. This is handy as a sanity check but is A TON of output. File itsxpress.stdout can be deleted after checking

`conda deactivate`

#### The following steps use three separate scripts to run the primary dada2 algorith, assign taxonomy, and then save ASV and taxonomy tables to hd along iwth some helpful graphs and tables. The workflow is split up to help avoid memory problems that can occur running the assignTaxonomy function on a lower memory system (with 8Gb of memory ~1K sequences can result in segfault and malloc errors). If running on a higher memory system the three separarte scripts can be run in one workflow with "dada2_full_workflow.r", or run serially as is done here.

Run dada2 algorith in R or run the script interactively
```
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/dada2.r
```
run taxonomy assignment
```
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/UNITE_taxonomic_classification.r
```
write tables to file and make plots
```
mkdir dada2_out
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/dada2_tables_to_file.r
```

The output includes asv table by counts in samples, asv taxonomy, asv rep seqs, asv seq lens by frequency table, sequences retained at various processing steps, as well as visualizations

### Combining run one `GARNAS_neonectria_barcoding_091819` and `GARNAS_neonectria_barcoding_012220`
Combine from ITS extraction step to run dada2 denoising on full dataset

```
mkdir GARNAS_neonectria_barcoding_runOneAndTwo_020320
mkdir GARNAS_neonectria_barcoding_runOneAndTwo_020320/R1_R2_switched/
mkdir GARNAS_neonectria_barcoding_runOneAndTwo_020320/R1_R2_switched/itsxpress

NEWDIR=GARNAS_neonectria_barcoding_runOneAndTwo_020320/R1_R2_switched/itsxpress
for i in GARNAS_neonectria_barcoding_091819/R1_R2_switched/itsxpress/*fastq.gz
do(
    FILE=${i##*/}
    cp $i ./$NEWDIR/runOne${FILE}
)
done

for i in GARNAS_neonectria_barcoding_012220/R1_R2_switched/itsxpress/*fastq.gz
do(
FILE=${i##*/}
cp $i ./$NEWDIR/runTwo${FILE}
)
done
```

Run dada2
```
cd GARNAS_neonectria_barcoding_runOneAndTwo_020320
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/dada2.r
Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/UNITE_taxonomic_classification.r
mkdir dada2_out
#Rscript ~/repo/neonectria_barcoding_012220/sequence_processing_and_dada2/dada2_tables_to_file.r

```

Ran comparison of dada2 ASV assignmnet in terms of combining samples across separate runs, and pooling parameters. These included concatenating sequence files from the same sample but different runs before ASV assigment, including all of the samples from the two separarte runs in a single `pool = T` call, and running the samples from two different runs through dada2 (with `pool = T`) separately. In the two latter approaches the sequences within an ASV were added from the same samples across the two runs. ASV richness, ASV overlap (between methods), and occurence of two Neonectria species was then compared. Concatenating the sample sequence files appear to be the most appropriate approach because 1) it creates the fewest number of ASVs that are unique to only that method, 2) it creates the fewest number of ASVs that occur in only a single sample, 3) it creates the fewest number of low abundance ASVs (i.e., it shortens rank abundance curve), and 4) it does not drop sequences compared to the other two methods. That is, the above suggests that we are culling sequence errors that would otherwise result in spurious richness.

### files_cat (i.e., files from the same sample concatenated and then pool=T in dada2) and files_sep (i.e. pool=T in dada2 alg. run across all files from both sequencing runs). These are performed after the intial sequence quality control and ITS extraction that is done on a per .fastq file basis.

```
Rscript ~/repo/neonectria_barcoding_012220/processing_methods_comparison/write_run1_run2_sample_pairs.r
#the bash script should be pointed at a folder of .fastq files that has been processed through sequence quality control and ITS extraction (itsxpress)
bash ~/repo/neonectria_barcoding_012220/processing_methods_comparison/cat_run1_run2_sample_pairs.sh
```
The next script was run on UNH HPC and performs several processing steps
```
sbatch ~/repo/neonectria_barcoding_012220/dada2_slurm.sh
```

