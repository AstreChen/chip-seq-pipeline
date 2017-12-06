#!/usr/bin/env python
import pandas as pd
import os,sys,subprocess
import argparse

parser = argparse.ArgumentParser(description='chipseq pipeline')
parser.add_argument('-c','--ctrl', nargs='*', help='input the ctrl sample name, e.g.: -c NAME, represents NAME.fastq.gz')
parser.add_argument('-t','--treat',nargs='*', help='input the case sample name, e.g.: -t NAME, represents NAME.fastq.gz')
parser.add_argument('-f','--factor',help='input the chip-seq antibody name.')
parser.add_argument('-r','--workdir',help='input the work dir, in this dir, a \'data/2.fastq/\' directory has been built for containing fastq files. please input the absolute path')

args = parser.parse_args()
treat = [x for x in args.treat]
ctrl = [x for x in args.ctrl]

#smpinfo = '/home/workdir/chen/proj/ighregion/geo40173/data/1.sample_info/p300_c-Myc_ctrl_info.csv'
refdir='/home/workdir/chen/genomes/mm10/mm10_bowtie2_idx/'

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

"""
#### ===== FASTQC ====== #####
for fname in os.listdir(rawdir) :
    cmd = 'fastqc '+ rawdir + '/' + fname +' -o '+fqcdir
    subprocess.call(cmd, shell=True)


#### ==== trim low quality reads === ####
for sname in SRRs:
    print sname
    cmd = 'sickle se -f '+ rawdir + sname+'.fastq.gz' + ' -t sanger -n -g -o ' + trimdir+sname+'.fastq.trimmed.gz ' \
        + '1>> '+logdir+'sickle.trimmed.logs 2>>'+logdir+'sickle.trimmed.errors'
    subprocess.call(cmd,shell=True)


### ==== FASTQC filter reads === ###
for fname in os.listdir(trimdir):
    cmd = 'fastqc '+ trimdir + '/' + fname +' -o '+fqcdir
    subprocess.call(cmd,shell=True)



### === Alignment === ###
for sname in SRRs:
    print sname
    cmd = 'bowtie2 -N 1 -x '+ refdir+'mm10 -U '+ trimdir+sname+'.fastq.trimmed.gz ' \
        + ' -S '+aligndir+sname+'.sam 1>> '+logdir+'bowtie2.aligns.logs 2>> '+ logdir+'bowtie2.aligns.errors ;' \
        + ' samtools view -S -b -o '+aligndir+sname+'.bam '+ aligndir+sname+'.sam ;' \
        + ' samtools sort '+aligndir+sname+'.bam '+aligndir+sname+'.sorted '
    subprocess.call(cmd,shell=True)
"""

### === call peak === ###
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
