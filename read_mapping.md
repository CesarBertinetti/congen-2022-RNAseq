# Read Mapping

## 1.	Index the reference using hisat2 build (Note: hisat2 executables must be in your PATH or specify the entire path to the executable)
```
cd reference
hisat2-build NW_025929828.fasta NW_025929828
cd ..
```

## 2.	Align the reads to the reference fasta using hisat2.
First it is necessary to make the directory that all data will be written to. note: this directory only needs to be created once.
```
mkdir alignments 

hisat2 -x reference/NW_025929828 -1 rawdata/CFA_S9.NW_025929828.R1.fastq.gz -2 rawdata/CFA_S9.NW_025929828.R2.fastq.gz -S alignments/CFA.sam
```
## 3.	Now that you have a sam file, you will want to convert the sam to a bam file, and check the mapping quality. Note: The output of samtools idxstats is tab-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.

Use samtools view to convert the SAM file into a BAM file and pipe to samtools sort. The following command requires samtools version 1.2 or higher.
```
cd alignments 
samtools view -bS CFA.sam | samtools sort - -o CFA.sort.bam
```

```
samtools index CFA.sort.bam
samtools idxstats CFA.sort.bam
cd ..
```
## 4.	Next use stringtie to generate a gtf for each sample for each gene in the reference annotation set on a per individual basis, this was an unstranded library preparation (it’s important to know how your data was generated), -e limits matches to the specified reference annotation file
```
stringtie alignments/CFA.sort.bam -G reference/GCF_023065955.1_UrsArc1.0_NW_025929828.gtf -e > CFA.gtf

mkdir ballgown
```
## 5.	Run for all samples. Let’s write a script that will do this for us! 
The shell script that works on my computer is saved as mapping_script.sh. 
Can use "touch mapping_script.sh" to create the blank file
```
input="sample_list.txt"
while IFS= read -r line
do
  echo "$line"
  hisat2 -x reference/NW_025929828 -1 rawdata/$line.NW_025929828.R1.fastq.gz -2 rawdata/$line.NW_025929828.R2.fastq.gz -S alignments/$line.sam
  samtools view -bS alignments/$line.sam | samtools sort - -o alignments/$line.sort.bam
  samtools index alignments/$line.sort.bam
  samtools idxstats alignments/$line.sort.bam > alignments/$line.stats
  stringtie alignments/$line.sort.bam -G reference/GCF_023065955.1_UrsArc1.0_NW_025929828.gtf -e -B -o ballgown/$line/$line.gtf
done < "$input"
```
To run the script, type
```
sh mapping_script.sh
```
At this point, you will have 12 bam files in separate folders, one for each individual 
## 6.	We want to generate counts data for each individual. There are many ways to do this, we will use prepDE.py, which is a python script provided with stringtie.  
```
prepDE_local.py -i ballgown
```
## 7.	We’ll now move to R for the remainder of our analyses. (will be provided as part2)
