# ChIP-Seq pipeline



## Motif Analysis
* Script: **motifAnalysis.py**
* Usage: 
 
 python motifAnalysis.py -p peakBEDfile -o OutputDir -n Name -g genome -find -annotate 

```
python motifAnalysis.py -p ctcf_narrowPeak -o ctcf_results -n ctcf -g hg19 -find -annotate
```

* Output:
 **name_motif1.bed** : bed file of motif in each peak region, including position, strand specific information.
 **name_motif1_annoPeak.txt**: bed file of each peak. Result format information can be found [here](http://homer.ucsd.edu/homer/ngs/peakMotifs.html) "Finding Instance of Specific Motifs"
 