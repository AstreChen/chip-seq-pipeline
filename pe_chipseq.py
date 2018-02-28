#!/usr/bin/env python
import pandas as pd
import os,sys,subprocess,fnmatch
import argparse
import logging

logging.basicConfig( level=logging.DEBUG)



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
def fastqc(SRRs, rawdir,fqcdir):
    logging.info("####Fastqc for quality control.")
    #for fname in os.listdir(rawdir) :
    for sname in SRRs:
        print sname
        files = fnmatch.filter(os.listdir(rawdir), sname+'*') 
        for f in files:
            cmd = 'fastqc '+ rawdir + '/' + f +' -o '+fqcdir
            subprocess.call(cmd, shell=True)
    logging.info("####Finish Fastqc.")


#### ==== trim low quality reads === ####
def trimReads(SRRs,rawdir,trimdir,logdir):
    logging.info("####Trimming low quality reads.")
    for sname in SRRs:
        print sname
        #cmd = 'sickle se -f '+ rawdir + sname+'.fastq.gz' + ' -t sanger -n -g -o ' + trimdir+sname+'.fastq.trimmed.gz ' \
        #    + '1>> '+logdir+'sickle.trimmed.logs 2>>'+logdir+'sickle.trimmed.errors'
        
        fqf1 = os.path.join(rawdir, '{}_R1.fastq.gz'.format(sname) )
        fqf2 = os.path.join(rawdir, '{}_R2.fastq.gz'.format(sname) )

        trmf1 = os.path.join(trimdir, '{}_R1.trimmed.fastq.gz'.format(sname) )
        trmf2 = os.path.join(trimdir, '{}_R2.trimmed.fastq.gz'.format(sname) )
        
        cmd = 'cutadapt -q 20 -m 20 -o {} -p {} {} {}'.format(trmf1, trmf2, fqf1, fqf2)
        #cmd = 'sickle pe -t sanger -s -g -f {} -r {} -o {} -p {} '.format( fqf1, fqf2, trmf1, trmf2) 
        #cmd2 = ' gzip {} {}'.format(trmf1, trmf2)
            #+ '1>> '+logdir+'sickle.trimmed.logs 2>>'+logdir+'sickle.trimmed.errors'
        print cmd
        #if not os.path.exists(trmf1):
        #    sys.exit("Failed in trimming reads.")

        subprocess.call(cmd  ,shell=True)

    logging.info("####Finished trimming.")


### ==== FASTQC filter reads === ###
def fqtrim(trimdir,fqcdir):
    logging.info("####Quality control after trimming reads.")
    for fname in os.listdir(trimdir):
        cmd = 'fastqc '+ trimdir + '/' + fname +' -o '+fqcdir
        subprocess.call(cmd,shell=True)

    logging.info("####Finished quality control.")



### === Alignment === ###
def alignReads(SRRs,refdir,args,trimdir,logdir,aligndir):
    logging.info("####Align reads...")

    refg = os.path.join(refdir, args.genome)

    for sname in SRRs:
        print sname
        trmf1 = os.path.join(trimdir, '{}_R1.trimmed.fastq.gz'.format(sname) )
        trmf2 = os.path.join(trimdir, '{}_R2.trimmed.fastq.gz'.format(sname) )
        samf = os.path.join(aligndir, sname+ '.sam' )
        bamf = os.path.join(aligndir, sname+'.bam' )
        sbamf = os.path.join(aligndir, sname+'.sorted.bam')

        aln_cmd = 'bowtie2 -N 1 -x {} -1 {} -2 {} -S {}'.format(refg, trmf1, trmf2, samf) + \
            ' 1>> '+logdir+'bowtie2.aligns.logs 2>> '+ logdir+'bowtie2.aligns.errors ;'

        s2bam_cmd = ' samtools view -S -b -o {} {}'.format(bamf, samf)
        bamsort_cmd = ' samtools sort {} -o {}'.format(bamf,sbamf)
        bamidx_cmd = ' samtools index {}'.format(sbamf)
        
        cmds = '\n'.join([aln_cmd, s2bam_cmd, bamsort_cmd, bamidx_cmd])
        print cmds
        subprocess.call(cmds,shell=True)
    logging.info("####Align reads complete.")

### === call peak === ###
def callPeak(args,aligndir,logdir,peakdir):
    logging.info("####Call ChIP peaks.")
    treat = [x for x in args.treat]
    ctrl = [x for x in args.ctrl]
    cmd = 'cd '+ aligndir + ' ; ' \
        + 'macs2 callpeak -t '+ '.sorted.bam '.join(treat) +'.sorted.bam ' \
        + ' -c '+ '.sorted.bam '.join(ctrl)+'.sorted.bam ' + ' --keep-dup 1 --nomodel -f BAMPE -g mm ' \
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
    logging.info("####Finished call peaks.")

if __name__=="__main__":
    fastqc(SRRs, rawdir,fqcdir)
    trimReads(SRRs,rawdir,trimdir,logdir)
    fqtrim(trimdir,fqcdir)
    alignReads(SRRs,refdir,args,trimdir,logdir,aligndir)
    callPeak(args,aligndir,logdir,peakdir)



