#!/bin/bash

#  cat_run1_run2_sample_pairs.sh
#  
#
#  Created by Eric Morrison on 2/17/20.
#  
#The files from run1 and run2 are being concatenated after itsxpress for processing through dada2 primary algorith. The intial quality control steps and itsxpress are performed on a read-wise basis, and so there is no need to process the reads from pairs of samples through the first part of the piepline as a sungle unit.


mkdir run1_run2_files_cat

while IFS= read -r line
do
    run1=$( echo "$line" | cut -f 2 )
    run2=$( echo "$line" | cut -f 3 )
    newLab=$( echo "$line" | cut -f 1 )

    run1R1=run1_run2_itsxpress_files_sep/${run1}_*_R1*.fastq.gz #Use negative match to 'P' to filter out sample repeats across plates; run 1 only
    run1R2=run1_run2_itsxpress_files_sep/${run1}_*_R2*.fastq.gz
    run2R1=run1_run2_itsxpress_files_sep/${run2}_*_R1*.fastq.gz
    run2R2=run1_run2_itsxpress_files_sep/${run2}_*_R2*.fastq.gz

    #echo "trying $run1R1"

    if [ -f $run1R1 ] && [ -f $run2R1 ]
    then
        #echo "success"
        #echo ""
        cat $run1R1 $run2R1 > run1_run2_files_cat/${newLab}_R1.fastq.gz
        cat $run1R2 $run2R2 > run1_run2_files_cat/${newLab}_R2.fastq.gz
    elif [ -f $run1R1 ] && [ ! -f $run2R1 ]
    then
        echo "$run1R1"
        cat $run1R1 > run1_run2_files_cat/${newLab}_R1.fastq.gz
        cat $run1R2 > run1_run2_files_cat/${newLab}_R2.fastq.gz
    elif [ ! -f $run1R1 ] && [ -f $run2R1 ]
    then
        echo "$run2R1"
        cat $run2R1 > run1_run2_files_cat/${newLab}_R1.fastq.gz
        cat $run2R2 > run1_run2_files_cat/${newLab}_R2.fastq.gz
    fi

done < sample_data/run1_run2_sample_pairs.txt
