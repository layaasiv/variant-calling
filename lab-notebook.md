# Lab Notebook - Variant Calling 

### Feb 22, 2025
Trying to roughly replicate this study: https://www.mdpi.com/2223-7747/9/4/439#app1-plants-09-00439 \
Found TAIR10 reference genome assembly on Ensembl. Downloaded via wget. \
Downloaded one of the study's DNAseq datasets using the SRR# provided and SRA toolkit. Here is the command used: \
```fasterq-dump --split-files SRR3340911```

For known sites VCF file, source was 1001 Genomes. 

To run gatk docker image: 

```docker run --rm -v /home/layaasiv/professional-dev/variant-calling/reference:/mnt broadinstitute/gatk gatk CreateSequenceDictionary -R /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta -O /mnt/Arabidopsis_thaliana.TAIR10.dna.toplevel.dict```


### Feb 23, 2025
Adapter seq for R1: (with index) GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGC; (without index) GATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
Adapter seq for R2: (index unknown) GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 

### Feb 24, 2025 

#### Running cutadapt

Command:
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o cutadapt/adapt_SRR3340911_1.fastq -p cutadapt/adapt_SRR3340911_2.fastq reads/SRR3340911_1.fastq reads/SRR3340911_2.fastq

Cutadapt output: 
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


#### Running trimmomatic

Command: trimmomatic -PE -threads 3 -phred33 cutadapt/adapt_SRR3340911_1.fastq cutadapt/adapt_SRR3340911_2.fastq \
trim_SRR3340911_1.fastq trim_unpaired_SRR3340911_1.fastq trim_SRR3340911_2.fastq trim_unpaired_SRR3340911_2.fastq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:85

Trimmomatic output:
Input Read Pairs: 17887298 
Both Surviving: 17556673 (98.15%) 
Forward Only Surviving: 194276 (1.09%) 
Reverse Only Surviving: 63219 (0.35%) 
Dropped: 73130 (0.41%)

Wanted to look at the distribution of read lengths after adaptor removal and trimming, but seems I'm not able to run the script to \
produce this data on local. I need to look for ways to get access to an HPC. 

## Feb 25, 2025 

#### BWA 
Created the reference index files for BWA. 
Command: bwa index ${ref}

Ran the BWA-MEM alignment, but forgot to add time command to survey how much resources it takes up, so canceled the run. It was taking \
a really long time to run though. Will start it again tomorrow. 

Command: bwa mem -t 4 ${ref} ${trimmed}/trim_SRR3340911_1.fastq ${trimmed}/trim_SRR3340911_2.fastq > ${aligned_reads}/SRR3340911_paired.sam

## Feb 27, 2025 
Ran MarkDuplicatesSpark successfully. 

When I tried running base recalibration, I ran into some errors: 

Error #1: 
```htsjdk.tribble.TribbleException$UnableToCreateCorrectIndexType: Unknown index type.  magic number: 0x49464356; type -1912514492```

This indicates that there is a problem with the known sites VCF file. I ended up downloading a new file for this from the TAIR10 ensembl page \
under the Variation section. I downloaded the VCF and CSI (index) files. 

Error #2: 
```java.lang.IllegalArgumentException: Something went wrong with sequence dictionary detection, check that reference has a valid sequence dictionary```

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
```A USER ERROR has occurred: Input /mnt/reference/arabidopsis_thaliana.vcf must support random access to enable queries by interval. If it's a file, please index it using the bundled tool IndexFeatureFile```

This indicates that there is an issue with the index file (.csi) for the known sites VCF file. To fix this, instead of using the .csi file, I created \
a new index file using the gatk IndexFeatureFile command: 

```
docker run --rm -v /home/layaasiv/professional-dev/variant-calling/reference:/mnt broadinstitute/gatk gatk IndexFeatureFile \
-I /mnt/arabidopsis_thaliana.vcf
```

Error #4:

```A USER ERROR has occurred: Error while trying to create index for /mnt/arabidopsis_thaliana.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 18: The VCF specification does not allow for whitespace in the INFO field. Offending field value was "The 1001 Genomes Project_2016;TSA=SNV"```

This indicates that there is an issue with the VCF file format, specifically that the known sites VCF file here contains spaces, which \
are invalid in the INFO field. Specifically, it's this entry: "The 1001 Genomes Project_2016;TSA=SNV". To fix this, I used ```sed``` to \
replace spaces with "_": 

```sed -i 's/ /_/g' arabidopsis_thaliana.vcf```

Error #5:
```htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 1687355: unparsable vcf record with allele W, for input source: /mnt/arabidopsis_thaliana.vcf```

This indicates there is some entry on line 1687355 of the known sites VCF file. Looking at that line specifcally with this command: \

```sed -n '1687355p' /mnt/arabidopsis_thaliana.vcf```

It indicates a W allele instead of one of the 4 regular bases, which is likely the reason for the error. I am not sure how to fix this \
entry, so I decided to remove the variant entirely with this command:

```sed -i '1687355d' /mnt/arabidopsis_thaliana.vcf```

Then, to ensure there were no other issues in the VCF file, I used the gatk ValidateVariants command:

```gatk ValidateVariants -V /mnt/arabidopsis_thaliana_fixed.vcf```

Finally, I re-ran the gatk BaseRecalibrator successfully.
