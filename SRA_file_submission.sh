#!/bin/sh

#  SRA_file_submission.sh
#  
#
#  Created by Eric Morrison on 2/17/21.
#  

mkdir ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs

for i in ~/GARNAS_neonectria_barcoding_091819/original_reads/*R1*
do(
    DIR=${i%/*}
    FILE=${i##*/}
    BEFFILE=${FILE%R1*}
    AFTFILE=${FILE##*R1}
    R2FILE=${BEFFILE}R2${AFTFILE}
    SAMPNAM=${FILE%%_*}

    cp $i ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs/runOne${SAMPNAM}_R1.fastq.gz
    cp $DIR/$R2FILE ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs/runOne${SAMPNAM}_R2.fastq.gz
)
done

for i in ~/GARNAS_neonectria_barcoding_012220/original_reads/*R1*
do(
DIR=${i%/*}
FILE=${i##*/}
BEFFILE=${FILE%R1*}
AFTFILE=${FILE##*R1}
R2FILE=${BEFFILE}R2${AFTFILE}
SAMPNAM=${FILE%%_*}

cp $i ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs/runTwo${SAMPNAM}_R1.fastq.gz
cp $DIR/$R2FILE ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs/runTwo${SAMPNAM}_R2.fastq.gz
)
done

for i in ~/GARNAS_neonectria_barcoding_files_cat_03242020/SRA_submission_seqs/*R1*
do(
DIR=${i%/*}
FILE=${i##*/}
R2FILE=${BEFFILE}R2${AFTFILE}
SAMPNAM=${FILE%%_*}

echo $SAMPNAM

)
done
