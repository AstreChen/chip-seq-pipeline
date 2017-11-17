#!/usr/bin/env bash
WORKDIR=/home/workdir/chen/proj/ighregion/embo2012/
LOG_DIR=$WORKDIR'/logs/'
DATA_DIR=$WORKDIR'/data/'
fastqc_dir=$DATA_DIR/fastqc
trimmed_dir=$DATA_DIR/trim

###================== Preprocess befor Alignment ====================#####
##### fastqc #####
#for fname in  $DATA_DIR/*.fastq.gz ; do
#    fastqc $fname -o $fastqc_dir
#done

### cut adapter and trim low quality reads ###
#for sname in `grep 'SRR' $WORKDIR/sample_SRR_info.txt | cut -f2` ; do
#    cutadapt -a 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTT' -g 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTT' -o $trimmed_dir/$sname'.trimmed.fastq.gz' $DATA_DIR/$sname'.fastq.gz' 1>> $LOG_DIR/cutadapt.trim.logs 2>> $LOG_DIR/cutadapt.trim.errors
#    sickle se -f $trimmed_dir/$sname'.trimmed.fastq.gz' -t sanger -n -g -o $trimmed_dir/$sname'.trimmed2.fastq.gz' 1>> $LOG_DIR/sickle.trimmed.logs
#done

### fastqc after filter reads ###
#for fname in $trimmed_dir/*.trimmed2.fastq.gz ; do
#    fastqc $fname -o $fastqc_dir
#done

### ==============  Alignment  ===================== ###
refdir='/home/workdir/chen/genomes/mm10/'
idxdir=$refdir'/mm10_bowtie2_idx/'
aligndir=$DATA_DIR/03.aligns
#for sname in `grep 'SRR' $WORKDIR/sample_SRR_info.txt | cut -f2` ; do
#    bowtie2 --local -N 0 -L 18 -x $idxdir/mm10 -U $trimmed_dir/$sname.trimmed2.fastq.gz -S $aligndir/$sname.sam 1>> $LOG_DIR/bowtie2.aligns.logs 2>> $LOG_DIR/bowtie2.aligns.errors
#    echo "###"$sname"\n"
#    samtools view -S -b -o $aligndir/$sname.bam $aligndir/$sname.sam
#    samtools sort $aligndir/$sname.bam $aligndir/$sname.sorted
#    samtools index $aligndir/$sname.sorted.bam $aligndir/$sname.idx
#done

peakdir=$DATA_DIR/04.peaks
cd $aligndir
macs2 callpeak -t \
    SRR499700.sorted.bam  SRR499702.sorted.bam  SRR499704.sorted.bam  SRR499706.sorted.bam \
    SRR499701.sorted.bam  SRR499703.sorted.bam  SRR499705.sorted.bam \
    -c SRR499707.sorted.bam SRR499708.sorted.bam \
    -f BAM -g mm --outdir $peakdir -n embo2012 -B -q 0.01 1>>$LOG_DIR/macs2.q01.logs  2>>$LOG_DIR/macs2.errors
    
