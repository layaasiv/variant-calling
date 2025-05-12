# Comparing Variant Calling Tools
Comparing different tools used in the variant calling pipeline. Replicating the work done in this paper: Comparison of Read Mapping and Variant Calling Tools for the Analysis of Plant NGS Data (Pucker et al. 2020).

## Aligners
Aligners being assessed in this project are: 
* BWA-MEM
* Bowtie2
* GEM3

### Comparing alignment metrics

![alt](/home/layaasiv/professional-dev/variant-calling/nf_revised/results_nf/metrics/plots/overall_alignment_metrics.png)

![alt](/home/layaasiv/professional-dev/variant-calling/nf_revised/results_nf/metrics/plots/alignment_error_metrics.png)

### Conclusions
BWA-MEM gives better alignment results than Bowtie2. 
Here are some possible reasons why Bowtie2 is performing poorly on this data: 

1) Bowtie2, by default, runs in "fast" mode, which trades off speed for sensitivity. BWA-MEM is sensitive by default.
2) Bowtie2 discards pairs more readily when insert size or orientation are not as expected.
3) Bowtie2 is also more strict about gaps and mismatches at default than BWA-MEM.
4) BWA-MEM also handles ambiguous/multi-mapping reads better, which is advantageous in Arabidopsis thaliana genome which has repetitive regions.

Overall, BWA-MEM gives better alignment results for genomic data from Arabidopsis thaliana. However, to improve the performance of Bowtie2, one method we can try is running it with the ```--very-sensitive``` parameter. Interestingly, the authors of the paper I'm replicating used only the default version of Bowtie2 and were able to get comparable mapping quality to BWA-MEM
