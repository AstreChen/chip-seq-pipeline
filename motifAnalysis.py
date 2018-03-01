#!/usr/bin/env python

import os,sys,subprocess,fnmatch
import argparse
import logging

logging.basicConfig( 
    level=logging.DEBUG, 
    format = '%(asctime)s %(filename)s : %(levelname)s  %(message)s',
    datefmt  = '%Y-%m-%d %A %H:%M:%S')

parser = argparse.ArgumentParser(description='Motif Analysis based on ChIP-seq')
parser.add_argument('-p','--peak', help='Input bed file like MACS2 output *.narrowPeak file.')
parser.add_argument('-o','--outdir', help='Output to this directory.')
parser.add_argument('-n','--name', help='Name of the factor.')
parser.add_argument('-g','--genome', help='Choose the genome, eg: hg19.')
parser.add_argument('-find', action='store_true', help='use -find to find all motifs. REQUIRED STEP!')
parser.add_argument('-annotate', action='store_true', help='use -annotate to annotate the best motif. Must use -find -annotate')

arg = parser.parse_args()

class MotifAnalysis:
    def __init__(self, peakfile, outdir, name, genome):
        self.peak = peakfile
        self.outdir = outdir
        self.name = name
        self.genome = genome
        if not os.path.exists(outdir):
            os.makedirs(outdir)


### Find motifs from narrowPeak ###

    # important parameter : -size 100
    def findMotifs(self):
        logging.info('Find motifs from region of input bed file .')
        cmd = 'findMotifsGenome.pl {} {} {} -size 100 -mask -len 10,12,14,16 -dumpFasta -bits -homer2'.format(self.peak, self.genome, self.outdir)
        subprocess.call(cmd, shell=True)
        logging.info('Find all motifs from region of input bed file. Output the result in the {}/homerResults.html. Please check the best motif ranked as first motif.')
        logging.info('!!!Finished finding all motif.')

### Annotate Peaks from best motif ###
    def annotateMotif(self):
        logging.info('Use the best motif to annotate all Peaks')
        cmd = 'annotatePeaks.pl {0} {1} -m {2}/homerResults/motif1.motif -mbed {2}/{3}_motif1.bed > {2}/{3}_motif1_annoPeak.txt'.format(self.peak, self.genome, self.outdir, self.name)
        subprocess.call(cmd, shell=True)
        logging.info('Annotate the position of peak region for the best predicted motif1. Output as {0}/{1}_motif1.bed and {0}/{1}_motif1_annoPeak.txt').format(self.outdir, self.name)
        logging.info('!!!Finished annotate best motif.')



if __name__ == '__main__':


    amotif = MotifAnalysis(peakfile=arg.peak, outdir=arg.outdir, name=arg.name, genome=arg.genome)
    if arg.find:
        amotif.findMotifs()
    if arg.annotate:
        amotif.annotateMotif()