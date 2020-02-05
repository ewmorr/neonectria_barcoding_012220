
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
Rscript ~/repo/neonectria_barcoding_012220/rm_primers_and_qual_filter.r
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
Rscript ~/repo/neonectria_barcoding_012220/dada2.r
```
run taxonomy assignment
```
Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification.r
```
write tables to file and make plots
```
mkdir dada2_out
Rscript ~/repo/neonectria_barcoding_012220/dada2_tables_to_file.r
```

The output includes asv table by counts in samples, asv taxonomy, asv rep seqs, asv seq lens by frequency table, sequences retained at various processing steps, as well as visualizations

