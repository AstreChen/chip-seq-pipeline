import pandas as pd
import numpy as np
import subprocess

outdir = '/home/workdir/chen/proj/ighregion/embo2012/data/'
df = pd.read_csv('../E-GEOD-38046.sdrf.txt',sep='\t')
cols = ['Source Name','Comment [Sample_source_name]',u'Comment [Sample_title]','Characteristics [organism part]','Comment [ENA_RUN]',u'Comment [ENA_RUN]',u'Comment [FASTQ_URI]',u'FactorValue [ANTIBODY]',u'FactorValue [EXTRACTION METHOD]', u'FactorValue [GENOTYPE]', u'FactorValue [ORGANISM PART]',u'FactorValue [STRAIN OR LINE]']

df_part = df.loc[:,cols]
dfchip=df_part.iloc[:,2].str.contains(pat='ChIP-seq')
dfchip = df_part.loc[dfchip,:]

## select pax5 ###
dfpax5 = dfchip.loc[dfchip.loc[:,'FactorValue [ANTIBODY]'].str.contains('Invitrogen'),:]
dfpax5.index = np.arange(dfpax5.shape[0])

for x in range(dfpax5.shape[0]):
    cmd = 'cd '+ outdir +';' + ' wget -c -t 0 '+ dfpax5.loc[x,'Comment [FASTQ_URI]']
    subprocess.call(cmd,shell=True)


