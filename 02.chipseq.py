#!/usr/bin/env python
import pandas as pd
import os,sys,subprocess
import argparse
import logging

logger = logging.getLogger('root')
logger.propagate = False


parser = argparse.ArgumentParser(description='chipseq pipeline')
parser.add_argument('-c','--ctrl', nargs='*', help='input the ctrl sample name, e.g.: -c NAME, represents NAME.fastq.gz')
parser.add_argument('-t','--treat',nargs='*', help='input the case sample name, e.g.: -t NAME, represents NAME.fastq.gz')
parser.add_argument('-f','--factor',help='input the chip-seq antibody name.')
parser.add_argument('-r','--workdir',help='input the work dir, in this dir, a \'data/2.fastq/\' directory has been built for containing fastq files. please input the absolute path')
parser.add_argument('-g','--genome',help='input the reference genome, "hg19" or "mm10", e.g. -g hg19. ')
args = parser.parse_args()


#smpinfo = '/home/workdir/chen/proj/ighregion/geo40173/data/1.sample_info/p300_c-Myc_ctrl_info.csv'
#refdir='/home/workdir/chen/genomes/mm10/mm10_bowtie2_idx/'
refdir='/home/ubuntu/genomes/bowtie2_indexes/'
treat = [x for x in args.treat]
ctrl = [x for x in args.ctrl]
SRRs = ctrl+treat
workdir = args.workdir +'/'
#workdir = '/home/workdir/chen/proj/ighregion/geo40173/'
logdir = workdir + 'logs/'
datadir = workdir + 'data/'
rawdir = datadir + '2.fastq/'
fqcdir = datadir + '/fastqc/'
trimdir = datadir + 'trim/'
aligndir = datadir + '3.align/'
peakdir = datadir + '4.peaks'

for folder in [logdir,rawdir,fqcdir,trimdir,aligndir,peakdir]:
    if not os.path.exists(folder):
        os.makedirs(folder)

#SRRs = pd.read_csv(smpinfo,sep='\t',header=None).iloc[:,0].tolist()

#### ===== FASTQC ====== #####
def fastqc(rawdir,fqcdir):
    logger.info("Fastqc for quality control.")
    for fname in os.listdir(rawdir) :
        cmd = 'fastqc '+ rawdir + '/' + fname +' -o '+fqcdir
        subprocess.call(cmd, shell=True)
    logger.info("Finish Fastqc.")


#### ==== trim low quality reads === ####
def trimReads(SRRs,rawdir,trimdir,logdir):
    logger.info("Trimming low quality reads.")
    for sname in SRRs:
        print sname
        cmd = 'sickle se -f '+ rawdir + sname+'.fastq.gz' + ' -t sanger -n -g -o ' + trimdir+sname+'.fastq.trimmed.gz ' \
            + '1>> '+logdir+'sickle.trimmed.logs 2>>'+logdir+'sickle.trimmed.errors'
        subprocess.call(cmd,shell=True)
    logger.info("Finished trimming.")


### ==== FASTQC filter reads === ###
def fqtrim(trimdir,fqcdir):
    for fname in os.listdir(trimdir):
        cmd = 'fastqc '+ trimdir + '/' + fname +' -o '+fqcdir
        subprocess.call(cmd,shell=True)



### === Alignment === ###
def alignReads(SRRs,refdir,args,trimdir,logdir,aligndir):
    logger.info("Align reads...")
    for sname in SRRs:
        print sname
        cmd = 'bowtie2 -N 1 -x '+ refdir+args.genome+' -U '+ trimdir+sname+'.fastq.trimmed.gz ' \
            + ' -S '+aligndir+sname+'.sam 1>> '+logdir+'bowtie2.aligns.logs 2>> '+ logdir+'bowtie2.aligns.errors ;' \
            + ' samtools view -S -b -o '+aligndir+sname+'.bam '+ aligndir+sname+'.sam ;' \
            + ' samtools sort '+aligndir+sname+'.bam '+aligndir+sname+'.sorted '
        subprocess.call(cmd,shell=True)
    logger.info("Align reads complete.")

### === call peak === ###
def callPeak(args,aligndir,logdir,peakdir):
    logger.info("Call ChIP peaks.")
    treat = [x for x in args.treat]
    ctrl = [x for x in args.ctrl]
    cmd = 'cd '+ aligndir + ' ; ' \
        + 'macs2 callpeak -t '+ '.sorted.bam '.join(treat) +'.sorted.bam ' \
        + ' -c '+ '.sorted.bam '.join(ctrl)+'.sorted.bam ' + ' --keep-dup 1 --nomodel -f BAM -g mm ' \
        + '--outdir '+peakdir+' -n '+args.factor +' -B -q 0.01 1>> '+logdir+'macs2_' + args.factor +'.logs 2>>'+logdir+'macs2_'+args.factor+'.errors '
    print cmd
    subprocess.call(cmd, shell=True)
    cmd = ' cd '+peakdir + '; LC_COLLATE=C sort -k1,1 -k2,2n '+ args.factor +'_treat_pileup.bdg > '+ args.factor +'_treat_pileup_sorted.bdg ;' \
    + '/home/workdir/chen/tools/bedGraphToBigWig '+ args.factor +'_treat_pileup_sorted.bdg /home/workdir/chen/genomes/mm10/chromInfo.txt.gz '+ args.factor +'.bw '
    print cmd
    subprocess.call(cmd, shell=True)

    cmd = ' cd '+peakdir + '; LC_COLLATE=C sort -k1,1 -k2,2n '+ args.factor +'_control_lambda.bdg > '+ args.factor +'_control_lambda_sorted.bdg ;' \
        + '/home/workdir/chen/tools/bedGraphToBigWig '+ args.factor +'_control_lambda_sorted.bdg /home/workdir/chen/genomes/mm10/chromInfo.txt.gz '+ args.factor +'_ctrl.bw '
    print cmd
    subprocess.call(cmd, shell=True)
    logger.info("Finished call peaks.")

if __name__=="__main__":
    fastqc(rawdir,fqcdir)
    trimReads(SRRs,rawdir,trimdir,logdir)
    fqtrim(trimdir,fqcdir)
    alignReads(SRRs,refdir,args,trimdir,logdir,aligndir)
    callPeak(args,aligndir,logdir,peakdir)



