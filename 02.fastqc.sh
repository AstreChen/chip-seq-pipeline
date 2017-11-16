#!/usr/bin/env bash
WORKDIR=/home/workdir/chen/proj/ighregion/embo2012/
LOG_DIR=$WORKDIR'/logs/'
DATA_DIR=$WORKDIR'/data/'
fastqc_dir=$DATA_DIR/fastqc
trimmed_dir=$DATA_DIR/trim

##### fastqc #####
#for fname in  $DATA_DIR/*.fastq.gz ; do
#    fastqc $fname -o $fastqc_dir
#done

#### trimmed low quality reads ###
for sname in `grep 'SRR' $WORKDIR/sample_SRR_info.txt | cut -f2` ; do
    fname=$sname'.fastq.gz'
    sickle se -f $DATA_DIR/$fname -t sanger -n -g -o $trimmed_dir/$sname'.fastq.trim.gz' 1>> $LOG_DIR/sickle.trim.logs
done
