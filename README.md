
mkdirs and rename R1 and R2 (this performs a "reverse complement" in essence, so that ITSxpress works properly)
```
cd GARNAS_neonectria_barcoding_091819


mkdir original_reads

for i in P*/*/*/*/reads/Sample*/*fastq.gz;do( mv $i ./original_reads/;)done

mkdir R1_R2_switched

for i in original_reads/*R1*; do(DIR=${i%/*}; FILE=${i##*/}; BEFFILE=${FILE%R1*}; AFTFILE=${FILE##*R1}; cp $DIR/${BEFFILE}R1${AFTFILE} ./R1_R2_switched/${BEFFILE}R2${AFTFILE}; cp $DIR/${BEFFILE}R2${AFTFILE} ./R1_R2_switched/${BEFFILE}R1${AFTFILE}); done

rm R1_R2_switched/*Empty*
```
Run the pre-processing script.
```
Rscript rm_primers_and_qual_filter.r --save --restore
```
Or open and R terminal and run the commands in `rm_primers_and_qual_filter.r` (This stuff should mostly run fine on a reasonably good laptop)

Next run itsxpress on the filtered reads. This should use the version of itsxpress that has been hacked to trim 3' read ends, other wise we end up with staggered read merges that have conserved 5.8S or LSU regions remaining (not good)

The `hmmsearch` step can take quite a while, so this is one place to use a cluster.

```
conda activate itsxpressHack

cd ~/GARNAS_neonectria_barcoding_091819/
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

    if [ -f $OUTDIR/$R2 ]
    then
        continue
    fi

    python ~/repo/itsxpress_mod/itsxpress/main.py \
    --fastq $INDIR/$R1 --fastq2 $INDIR/$R2 \
    --outfile $OUTDIR/$R1 --outfile2 $OUTDIR/$R2 \
    --region ITS2 --taxa 'Fungi' --cluster_id 1 \
    --threads 4 \
    --log itsxpress.log >> itsxpress.stdout
)
done
```
The stdout is sent to file bc I have the program set to output the trimmed read len of each fwd and reverse read. This is handy as a sanity check but is A TON of output. File itsxpress.stdout can be deleted after checking
conda deactivate

After itsxpress 237 of the original 280 read file pairs is remaining, (239 were input to itsxpress).

Run dada2 algorith in R
```
Rscript dada2.r --save --restore
```
or run the script interactively

The output includes asv table by counts in samples, asv taxonomy, asv rep seqs, asv seq lens by frequency table, sequences retained at various processing steps, as well as visualizations

