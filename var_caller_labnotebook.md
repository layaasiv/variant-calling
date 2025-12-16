# Lab Notebook - Variant Calling 

## Feb 22, 2025
Trying to roughly replicate this study: https://www.mdpi.com/2223-7747/9/4/439#app1-plants-09-00439 \
Found TAIR10 reference genome assembly on Ensembl. Downloaded via wget. \
Downloaded one of the study's DNAseq datasets using the SRR# provided and SRA toolkit. Here is the command used: \

```fasterq-dump --split-files SRR3340911```

For known sites VCF file, source was 1001 Genomes. 

To run gatk docker image: 

```
docker run --rm -v /home/layaasiv/professional-dev/variant-calling/reference:/mnt broadinstitute/gatk gatk CreateSequenceDictionary -R /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta -O /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.dict
```


## Feb 23, 2025
Adapter seq for R1: (with index) GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC; (without index) GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
Adapter seq for R2: (index unknown) GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 

## Feb 24, 2025 

#### Running cutadapt

Command:

```
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o cutadapt/adapt_SRR3340911_1.fastq -p cutadapt/adapt_SRR3340911_2.fastq reads/SRR3340911_1.fastq reads/SRR3340911_2.fastq
```

Cutadapt output: 

```
=== Summary ===

Total read pairs processed:         17,887,298
  Read 1 with adapter:                 281,320 (1.6%)
  Read 2 with adapter:                 236,438 (1.3%)
Pairs written (passing filters):    17,887,298 (100.0%)

Total basepairs processed: 8,943,649,000 bp
  Read 1: 4,471,824,500 bp
  Read 2: 4,471,824,500 bp
Total written (filtered):  8,930,212,125 bp (99.8%)
  Read 1: 4,459,850,413 bp
  Read 2: 4,470,361,712 bp
```

#### Running trimmomatic

Command: 

```
trimmomatic PE -threads 3 -phred33 cutadapt/adapt_SRR3340911_1.fastq cutadapt/adapt_SRR3340911_2.fastq \
trim_SRR3340911_1.fastq trim_unpaired_SRR3340911_1.fastq trim_SRR3340911_2.fastq trim_unpaired_SRR3340911_2.fastq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:85
```

Trimmomatic output:

```
Input Read Pairs: 17887298
Both Surviving: 17556673 (98.15%)
Forward Only Surviving: 194276 (1.09%) 
Reverse Only Surviving: 63219 (0.35%) 
Dropped: 73130 (0.41%)
```

Wanted to look at the distribution of read lengths after adaptor removal and trimming, but seems I'm not able to run the script to \
produce this data on local. I need to look for ways to get access to an HPC. 

## Feb 25, 2025 

#### BWA 
Created the reference index files for BWA. 

Command: 

```
bwa index ${ref}
```

Ran the BWA-MEM alignment, but forgot to add time command to survey how much resources it takes up, so canceled the run. It was taking \
a really long time to run though. Will start it again tomorrow. 

Command: 

```
bwa mem -t 4 ${ref} ${trimmed}/trim_SRR3340911_1.fastq ${trimmed}/trim_SRR3340911_2.fastq > ${aligned_reads}/SRR3340911_paired.sam
```

## Feb 27, 2025 
Ran MarkDuplicatesSpark successfully. 

When I tried running base recalibration, I ran into some errors: 

Error #1: 

```
htsjdk.tribble.TribbleException$UnableToCreateCorrectIndexType: Unknown index type.  magic number: 0x49464356; type -1912514492
```

This indicates that there is a problem with the known sites VCF file. I ended up downloading a new file for this from the TAIR10 ensembl page \
under the Variation section. I downloaded the VCF and CSI (index) files. 

Error #2: 

```
java.lang.IllegalArgumentException: Something went wrong with sequence dictionary detection, check that reference has a valid sequence dictionary
```

This indicates that there is an issue with the sequence dictionary file that we previously created for the reference genome. A properly formatted \
.dict file should look something like this: 

```
@HD    VN:1.0  SO:coordinate
@SQ    SN:chr1  LN:248956422
@SQ    SN:chr2  LN:242193529
```

But the file I created had only 1 line: 

```
@HD     VN:1.6
```

So, I re-ran the gatk CreateSequenceDictionary command:

```
docker run --rm -v /home/layaasiv/professional-dev/variant-calling/reference:/mnt broadinstitute/gatk gatk CreateSequenceDictionary \
-R /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
-O /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.dict
```

Error #3:

```
A USER ERROR has occurred: Input /mnt/reference/arabidopsis_thaliana.vcf must support random access to enable queries by interval. If it's a file, please index it using the bundled tool IndexFeatureFile
```

This indicates that there is an issue with the index file (.csi) for the known sites VCF file. To fix this, instead of using the .csi file, I created \
a new index file using the gatk IndexFeatureFile command: 

```
docker run --rm -v /home/layaasiv/professional-dev/variant-calling/reference:/mnt broadinstitute/gatk gatk IndexFeatureFile \
-I /mnt/arabidopsis_thaliana.vcf
```

Error #4:

```
A USER ERROR has occurred: Error while trying to create index for /mnt/arabidopsis_thaliana.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 18: The VCF specification does not allow for whitespace in the INFO field. Offending field value was "The 1001 Genomes Project_2016;TSA=SNV"
```

This indicates that there is an issue with the VCF file format, specifically that the known sites VCF file here contains spaces, which \
are invalid in the INFO field. Specifically, it's this entry: "The 1001 Genomes Project_2016;TSA=SNV". To fix this, I used sed to \
replace spaces with "_": 

```
sed -i 's/ /_/g' arabidopsis_thaliana.vcf
```

Error #5:

```
htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 1687355: unparsable vcf record with allele W, for input source: /mnt/arabidopsis_thaliana.vcf
```

This indicates there is some entry on line 1687355 of the known sites VCF file. Looking at that line specifcally with this command: \

```
sed -n '1687355p' /mnt/arabidopsis_thaliana.vcf
```

It indicates a W allele instead of one of the 4 regular bases, which is likely the reason for the error. I am not sure how to fix this \
entry, so I decided to remove the variant entirely with this command:

```
sed -i '1687355d' /mnt/arabidopsis_thaliana.vcf
```

Then, to ensure there were no other issues in the VCF file, I used the gatk ValidateVariants command:

```
gatk ValidateVariants -V /mnt/arabidopsis_thaliana_fixed.vcf
```

Finally, I re-ran the gatk BaseRecalibrator successfully.


## Feb 28, 2025

Ran gatk ApplyBQSR successfully:

```
Percent of CPU this job got: 0%
Elapsed (wall clock) time (h:mm:ss or m:ss): 26:30.54
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 37456
```

Ran gatk CollectAlignmentMetrics successfully:

```
Percent of CPU this job got: 0%
Elapsed (wall clock) time (h:mm:ss or m:ss): 4:49.63
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 34452
```

Ran gatk CollectInsertSizeMetrics successfully:

```
Percent of CPU this job got: 0%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:52.95
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 33608
```

Ran gatk HaplotypeCaller successfully:

```
Percent of CPU this job got: 0%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:50:36
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 35080
```


## Mar 1, 2025


Ran Bowtie 2 aligner. 

Output:

```
17556673 reads; of these:
  17556673 (100.00%) were paired; of these:
    15185087 (86.49%) aligned concordantly 0 times
    1851602 (10.55%) aligned concordantly exactly 1 time
    519984 (2.96%) aligned concordantly >1 times
    ----
    15185087 pairs aligned concordantly 0 times; of these:
      81882 (0.54%) aligned discordantly 1 time
    ----
    15103205 pairs aligned 0 times concordantly or discordantly; of these:
      30206410 mates make up the pairs; of these:
        21932819 (72.61%) aligned 0 times
        6010774 (19.90%) aligned exactly 1 time
        2262817 (7.49%) aligned >1 times
37.54% overall alignment rate
```


```
Percent of CPU this job got: 19%
Elapsed (wall clock) time (h:mm:ss or m:ss): 13:48:48
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 174204
```

For issues with accessing docker: ```sudo usermod -aG docker $USER```


## Mar 29, 2025

Trying to convert the variant calling pipeline into a nextflow script. 
Make sure you're in the nextflow env before running nf pipelines because that's where Java 11 is installed.

Was previously running into issues with CreateSeqDict and FastQC processes. So, I isolated fastqc into a separate nf script to troubleshoot and fixed it. It now runs successfully. 

There were a couple of issues I noticed:

1) 1 param variable was defined twice 
2) The directory you want for the results to be stored and accessed outside the working dir should be specified in publishDir in the process. However, in the 'oputput' part of the process definition, you just mention dir which will be in the working dir. In the script, before you run the fastqc, mkdir for the dir mentioned in 'ouput'. Overall, do not mix the definition of these 2 components up within the process definition.

Next, to troubleshoot CreateSeqDict process. 

## Apr 4, 2025

Got the nextflow pipeline working and better understanding how this works. 

Trying to make same progress with the bowtie2 pipeline. Running MarkDuplicatesSpark and BQSR, but got error that the sam file produced by bowtie2 does not contain RG lines. Therefore, had to install picard package: 

```conda install bioconda::picard```

And then run this function:

```
picard AddOrReplaceReadGroups \
I=bowtie_aligned_SRR3340911.sam \
O=bowtie_aligned_SRR3340911_rg.sam \
RGID=SRR3340911 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=seq \
RGSM=SRR3340911
```
Elapsed time: 3.72 minutes.
Runtime.totalMemory()=715128832

RGID: can be set to anything that makes sense \
RGLB: library -- library prep kit used \
RGPL: platform \
RGPU: platform unit -- i.e., flowcell or sequencing run \
RGSM: sample name

Initially, tried to run with just the ID, PL and SM info (as in BWA), but the other fields were also required in this case. So, inserted some arbitrary terms. 


## May 1, 2025 
Successfully running the nextflow script now. 

I have combined both BWA and bowtie2 processes into the same pipeline. 

Comparing the aligners using these metrics:

1) samtools flagstat
2) gatk collect alignment metrics 
3) gatk collect insert size metrics 

Found that BWA in general gives better alignment results for this data. Here is a summary of the results extracted from the samtools flagstat output:

| Metric | BWA | Bowtie2 | Notes |
| ------ | --- | ------- | ----- |
| Total reads | 52.5M | 35.1M | Seems like Bowtie2 has been run on a subset of reads, even though this is not the case |
| Mapped (%) | 99.10% | 37.54% | Massive difference here — very low for Bowtie2 |
| Properly paired (%) | 88.48% | 13.51%	| Huge quality gap |
| Duplicates | 20.4M | 12.3M | Expected — scales with mapped reads |
| Singletons | 0.15% | 23.31% | Indicates many reads in Bowtie2 had unmapped mates |
| Supplementary | 17.3M | 0 | BWA produced lots of supplementary alignments (e.g. split reads, chimeras) — Bowtie2 didn’t report any |


Since both BWA-MEM and Bowtie2 used the same reference files, these likely are not an issue. \
Here are some possible reasons why Bowtie2 is performing poorly on this data: 

1) Bowtie2, by default, runs in "fast" mode, which trades off speed for sensitivity. BWA-MEM is sensitive by default.
2) Bowtie2 discards pairs more readily when insert size or orientation are not as expected.
3) Bowtie2 is also more strict about gaps and mismatches at default than BWA-MEM.
4) BWA-MEM also handles ambiguous/multi-mapping reads better, which is advantageous in Arabidopsis thaliana genome which has repetitive regions.

Overall, BWA-MEM gives better alignment results for genomic data from Arabidopsis thaliana. However, to improve the performance of Bowtie2, one method we can try is running it with the ```--very-sensitive``` parameter. Interestingly, the authors of the paper I'm replicating used only the default version of Bowtie2 and were able to get comparable mapping quality to BWA-MEM. 

## May 9, 2025
Also added GEM3 to the analysis. I modularized the nextflow script so it is easier to add new processes and parallelize the same process on multiple files. 

Looking at the alignment metrics now. Compiling them all into a multiqc report, and isolated plots of interest to present in the README.md. 

Used this tool to get insert size metrics:

```
picard CollectInsertSizeMetrics \
I=your_aligner_output.bam \
O=insert_metrics.txt \
H=insert_histogram.pdf \
M=0.5
```

Made an insert size histogram. Found that BWA-MEM has a bimodal distribution for the insert sizes, which is abnormal. We expect to see a unimodal, normal-looking curve.
Looking at the values from picard CollectInsertSizeMetrics, we see a max insert size of 29,008,098 (BWA) and 16,778,582 (bowtie2), and a min insert size of 2 (BWA) and 75 (bowtie2). These high insert size values are likely the result of mis-paired reads aligning to repetitive regions. The very small insert sizes could be the result of adapter contamination during sequencing, incomplete adapter ligation during library prep or other low complexity artifacts. Since variant callers usually expect normally distributed insert size metrics, these outliers should be filtered beforehand.

This will be accomplished using bamtools:

```
bamtools filter \
-in markeddup.bam \
-insertSize ">=100" \
-insertSize "<=2000" \
-out filtered.bam
```

```
samtools view -h SRR3340911_bwa_markeddup.bam | \
awk '($1 ~ /^@/) || ($9 >= 100 && $9 <= 2000) || ($9 <= -100 && $9 >= -2000)' | \
samtools view -b -o SRR3340911_bwa_insert_filtered_abs.bam -
```

## May 10, 2025
I filtered for insert sizes, and re-ran the metrics. Working on the plots now. 

## May 22, 2025
I completed analysis of aligners. BWA-MEM was the best, so decided to move forward with the BAM output from that for variant callers comparison. I chose to stick to 2 variant callers: GATK HaplotypeCaller and FreeBayes. 

Ran the var callers and did evaluation using Illumina hap.py. This required a ground truth vcf file. In the original paper, they linked a gold standard Nd1 vs TAIR10 variants vcf file that was used; I used the same one. The raw file did not have a header so I had to add that information in. Also, there were some bases that were unrecognized, so I filtered the vcf to retain only those entries that contained ACTG. 

The paper also mentioned that they used a Nd1-specific reference genome for the evaluation part of the experiment. I downloaded this reference, but was unable to use it in hap.py because the contig lengths did not match between the query and reference (since the query vcf were generated using the TAIR10 reference genome). So, I used TAIR10 as the reference in the evaluation as well. This also made the most sense to me since it maintains consistency throughout the experiment; it is not clear to me why they switched the reference. 

They also seem to have used their own evaluation pipeline (not sure if it included hap.py). 

The results were ultimately very similar between the variant callers, however HaplotypeCaller has better precision and F1 score in this experiment. These results, of course, are specific to this one sample. 

Future directions would include having more samples (replicates) and also exploring the performance of these tools on sequencing data from different organisms. Also, would test more tools; was limited here by compute resources.  