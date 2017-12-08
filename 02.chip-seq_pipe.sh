#!/usr/bin/env bash
WORKDIR=/home/workdir/chen/proj/ighregion/geo40173
LOG_DIR=$WORKDIR'/logs/'
DATA_DIR=$WORKDIR'/data/'
raw_dir=$DATA_DIR/2.fastq
fastqc_dir=$DATA_DIR/fastqc
trimmed_dir=$DATA_DIR/trim
aligndir=$DATA_DIR/3.aligns
peakdir=$DATA_DIR/4.peaks

for folder in ($LOG_DIR,$fastqc_dir,$trimmed_dir,$aligndir,$peakdir)
    if [! -x "$folder"]; then
        mkdir "$folder"
    fi

SRRs=`cut -f1 /home/workdir/chen/proj/ighregion/geo40173/data/1.sample_info/p300_c-Myc_ctrl_info.csv`
###================== Preprocess befor Alignment ====================#####
##### fastqc #####
for fname in  $raw_dir/*.fastq.gz ; do
    fastqc $fname -o $fastqc_dir
done

### cut adapter and trim low quality reads ###
#for sname in $SRRs; do
#    cutadapt -a 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTT' -g 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTT' -o $trimmed_dir/$sname'.trimmed.fastq.gz' $DATA_DIR/$sname'.fastq.gz' 1>> $LOG_DIR/cutadapt.trim.logs 2>> $LOG_DIR/cutadapt.trim.errors
#    sickle se -f $raw_dir/$sname'.trimmed.fastq.gz' -t sanger -n -g -o $trimmed_dir/$sname'.trimmed2.fastq.gz' 1>> $LOG_DIR/sickle.trimmed.logs
#done

### fastqc after filter reads ###
#for fname in $trimmed_dir/*.trimmed2.fastq.gz ; do
#    fastqc $fname -o $fastqc_dir
#done

### ==============  Alignment  ===================== ###
refdir='/home/workdir/chen/genomes/mm10/'
idxdir=$refdir'/mm10_bowtie2_idx/'
#for sname in `grep 'SRR' $WORKDIR/sample_SRR_info.txt | cut -f2` ; do
#    bowtie2 --local -N 0 -L 18 -x $idxdir/mm10 -U $trimmed_dir/$sname.trimmed2.fastq.gz -S $aligndir/$sname.sam 1>> $LOG_DIR/bowtie2.aligns.logs 2>> $LOG_DIR/bowtie2.aligns.errors
#    echo "###"$sname"\n"
#    samtools view -S -b -o $aligndir/$sname.bam $aligndir/$sname.sam
#    samtools sort $aligndir/$sname.bam $aligndir/$sname.sorted
#    samtools index $aligndir/$sname.sorted.bam $aligndir/$sname.idx
#done

#cd $aligndir
#macs2 callpeak -t \
#    -c SRR499707.sorted.bam SRR499708.sorted.bam \
#    --SPMR --keep-dup 1 --nomodel \
#    -f BAM -g mm --outdir $peakdir -n pax5 -B -q 0.01 1>>$LOG_DIR/macs2.1120.logs  2>>$LOG_DIR/macs2.1120.errors

#macs2 bdgcmp -t embo2012_treat_pileup.bdg -c embo2012_control_lambda.bdg -o pax5_FE.bdg -m FE
#LC_COLLATE=C sort -k1,1 -k2,2n pax5_FE.bdg > pax5_FE_sorted.bdg
#/home/workdir/chen/tools/bedGraphToBigWig pax5_FE_sorted.bdg /home/workdir/chen/genomes/mm10/chromInfo.txt.gz pax5_FE.bwig

