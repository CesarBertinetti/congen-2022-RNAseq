# congen-2022
Code and descriptions for ConGen 2022

RNAseq analysis methods
Prerequisites:
* Python 2.7
* R
* EdgeR (installation described below within R)
* HISAT2 v2.2.0 https://ccb.jhu.edu/software/hisat2/index.shtml
* StringTie v2.1.4 https://ccb.jhu.edu/software/stringtie/index.shtml#install
* prepDE.py https://ccb.jhu.edu/software/stringtie/dl/prepDE.py
* Samtools v1.10 http://www.htslib.org/download/

Worksheet: 
We have RNAseq data from Jansen et al https://www.nature.com/articles/s42003-019-0574-4. The reads have already been adapter trimmed and quality trimmed. The necessary reference files (gtf and fasta) have been downloaded and are available in the data folder. The files were downloaded from: 

https://www.ncbi.nlm.nih.gov/genome/23993?genome_assembly_id=1832096

For this tutorial, I selected only chromosome and annotations from NW_025929828. 

A bit more about gtf format:
http://www.ensembl.org/info/website/upload/gff.html

