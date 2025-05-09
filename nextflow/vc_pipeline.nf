/*
 * Pipeline parameters
 */

// Primary input
params.ref_fasta = "/home/layaasiv/professional-dev/variant-calling/nextflow/ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
params.read1 = "/home/layaasiv/professional-dev/variant-calling/reads/SRR3340911_1.fastq"
params.read2 = "/home/layaasiv/professional-dev/variant-calling/reads/SRR3340911_2.fastq"
params.refvcf = "/home/layaasiv/professional-dev/variant-calling/reference/arabidopsis_thaliana.vcf"

// Output directory
params.outdir = "/home/layaasiv/professional-dev/variant-calling/nextflow/results_nf"

/*
 * Generate fasta index file
 */

process samtools_index {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir params.outdir, mode: "symlink"

    input:
        path input_fasta

    output:
        path "${input_fasta}.fai", emit: "ref_fasta_fai"

    script:
    """
    samtools faidx "${input_fasta}"
    """

}

/*
 * Create sequence dictionary for the reference genome fasta file
 */

process gatk_createseqdict {
    container "broadinstitute/gatk"

    publishDir params.outdir, mode: "symlink"

    input:
        path input_fasta
    
    output:
        path "${input_fasta.baseName}.dict", emit: "ref_fasta_dict"

    script:
    """
    gatk CreateSequenceDictionary \
    -R ${input_fasta} \
    -O ${input_fasta.baseName}.dict
   
    """
}

/*
 * Create index file for the reference VCF feature file
 */
process gatk_indexvcf {
    container "broadinstitute/gatk"

    publishDir params.outdir, mode: "symlink"

    input:
        path ref_vcf
    
    output:
        path "${ref_vcf}.idx", emit: "ref_vcf_idx"
    
    script:
    """
    gatk IndexFeatureFile \
    -I ${ref_vcf}
    """
}

/*
 * Run FastQC on the test reads 
 */

process fastqc {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir params.outdir, mode: "symlink"

    input:
        path input_read1
        path input_read2

    output:
        path "fastqc_reports/*"
        path "fastqc_reports/*"

    script:
    """
    # Create the output directory if it doesn"t exist
    mkdir -p fastqc_reports
    
    # run fastqc
    fastqc ${input_read1} --outdir fastqc_reports/

    fastqc ${input_read2} --outdir fastqc_reports/
    """
}

/*
 * Remove adapter sequences from reads using cutadapt
 */

process cutadapt {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

    input:
        path input_read1
        path input_read2

    output:
        path "adapt_${input_read1}", emit: "adapt_read1"
        path "adapt_${input_read2}", emit: "adapt_read2"
    
    script:
    """
    cutadapt \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o adapt_${input_read1} \
    -p adapt_${input_read2} \
    ${input_read1} ${input_read2}
    """
}

/*
 * Trim reads using trimmomatic
 */

process trimmomatic {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

    input:
        path cutadapt_output_read1
        path cutadapt_output_read2

    output:
        path "trim_read1.fastq", emit: "trimmed_read1"
        path "trim_unpaired_read1.fastq"
        path "trim_read2.fastq", emit: "trimmed_read2"
        path "trim_unpaired_read2.fastq"
    
    script:
    """
    trimmomatic PE -threads 3 -phred33 \
    ${cutadapt_output_read1} ${cutadapt_output_read2} \
    trim_read1.fastq trim_unpaired_read1.fastq \
    trim_read2.fastq trim_unpaired_read2.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:85
    """
}


/*
 * Aligning the reads to reference genome using BWA
 */

process BWA {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/bwa", mode: "symlink"

    input:
        path ref_fasta
        path trimmed_read1
        path trimmed_read2
        val base_fn
    
    output:
        path "${ref_fasta}.amb"
        path "${ref_fasta}.ann"
        path "${ref_fasta}.bwt"
        path "${ref_fasta}.pac"
        path "${ref_fasta}.sa"
        tuple path("*.sam"), val(base_fn) into bwa_out_ch
    
    script:
    """
    echo "Creating index files for reference for BWA usage"
    bwa index ${ref_fasta}

    echo "Running BWA alignment"
    bwa mem -t 4 -R "@RG\\tID:SRR3340911\\tPL:ILLUMINA\\tSM:SRR3340911" \
    ${ref_fasta} \
    ${trimmed_read1} ${trimmed_read2} \
    > ${base_fn}_bwa.sam
    """
}

/*
 Run Bowtie2 aligner
 */

process bowtie {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/bowtie", mode: "symlink"

    input:
        path ref_fasta
        val idx_filename_base
        path trimmed_read1
        path trimmed_read2
        val base_fn

    output:
        path "${idx_filename_base}.1.bt2"
        path "${idx_filename_base}.2.bt2"
        path "${idx_filename_base}.3.bt2"
        path "${idx_filename_base}.4.bt2"
        path "${idx_filename_base}.rev.1.bt2"
        path "${idx_filename_base}.rev.2.bt2"
        path "${base_fn}_bowtie_raw.sam", emit: "bt_sam"

    script:
    """
    echo "Building reference index files for bowtie2 aligner"

    bowtie2-build \
    ${ref_fasta} \
    ${idx_filename_base}

    echo "Running bowtie2 aligner"

    bowtie2 \
    -x ${idx_filename_base} \
    -1 ${trimmed_read1} \
    -2 ${trimmed_read2} > \
    ${base_fn}_bowite_raw.sam
    """
}

/*
 Adding read group information to the Bowtie .sam file
 */

process add_rg {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/bowtie", mode: "symlink"

    input:
        path bowtie_aligned_sam
        val RGID
        val RGLB
        val RGPL
        val RGPU
        val RGSM
        val base_fn
    
    output:
        tuple path("*.sam"), val(base_fn) into bowtie_out_ch

    script:
    """
    picard AddOrReplaceReadGroups \
    I=${bowtie_aligned_sam} \
    O=${base_fn}_bowtie.sam \
    RGID=${RGID} \
    RGLB=${RGLB} \
    RGPL=${RGPL} \
    RGPU=${RGPU} \
    RGSM=${RGSM}
    """
}

/*
Sort SAM files
*/

process st_sort_index {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/metrics", mode: "symlink"

    input:
        path raw_aligned_sam
        val out_bam_fn

    output:
        path "${out_bam_fn}_sorted.bam"
        path "${out_bam_fn}_sorted.bam.bai"

    script:
    """
    samtools -Sb ${raw_aligned_sam} | samtools sort -o ${out_bam_fn}_sorted.bam

    samtools index ${out_bam_fn}_sorted.bam
    """
}

/*
 * Mark duplicates & sort using GATK
 */

process markdupsort {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}", mode: "symlink"

    input:
        path bwa_aligned_sam
        val bwa_output_filename
        path bt_aligned_sam
        val bt_output_filename

    output:
        path "bwa/${bwa_output_filename}_sorted_dedup_reads.bam", emit: "bwa_sort_dedup_bam"
        path "bowtie/${bt_output_filename}_sorted_dedup_reads.bam", emit: "bt_sort_dedup_bam"

    script:
    """
    mkdir -p bowtie

    echo "Marking duplicates for BWA aligned reads"

    gatk MarkDuplicatesSpark \
    -I ${bwa_aligned_sam} \
    -O bwa/${bwa_output_filename}_sorted_dedup_reads.bam

    echo "Marking duplicates for Bowtie aligned reads"

    gatk MarkDuplicatesSpark \
    -I ${bt_aligned_sam} \
    -O bowtie/${bt_output_filename}_sorted_dedup_reads.bam
    """
}

/*
 * Base quality score recalibration
 */

process bqsr {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}", mode: "symlink"

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path ref_vcf
        path ref_vcf_idx
        path bwa_sorted_dedup_reads_bam
        val bwa_output_recal_table
        path bt_sorted_dedup_reads_bam
        val bt_output_recal_table

    output:
        path "bwa/*.table", emit: "bwa_bqsr_recal_table"
        path "bwa/*_bqsr.bam", emit: "bwa_bqsr_bam"
        path "bowtie/*.table", emit: "bt_bqsr_recal_table"
        path "bowtie/*_bqsr.bam", emit: "bt_bqsr_bam"

    script:
    """
    mkdir -p bwa bowtie

    echo "Building the BQSR model: BWA"

    gatk BaseRecalibrator \
    -I ${bwa_sorted_dedup_reads_bam} \
    -R ${ref_fasta} \
    --known-sites ${ref_vcf} \
    -O bwa/${bwa_output_recal_table}.table

    echo "Applying the BQSR model: BWA"

    gatk ApplyBQSR \
    -I ${bwa_sorted_dedup_reads_bam} \
    -R ${ref_fasta} \
    --bqsr-recal-file bwa/${bwa_output_recal_table}.table \
    -O bwa/${bwa_sorted_dedup_reads_bam.baseName}_bqsr.bam

    echo "Building the BQSR model: Bowtie"

    gatk BaseRecalibrator \
    -I ${bt_sorted_dedup_reads_bam} \
    -R ${ref_fasta} \
    --known-sites ${ref_vcf} \
    -O bowtie/${bt_output_recal_table}.table

    echo "Applying the BQSR model: Bowtie"

    gatk ApplyBQSR \
    -I ${bt_sorted_dedup_reads_bam} \
    -R ${ref_fasta} \
    --bqsr-recal-file bowtie/${bt_output_recal_table}.table \
    -O bowtie/${bt_sorted_dedup_reads_bam.baseName}_bqsr.bam
    """
}

/*
 Collect alignment and insert size metrics
 */

process gatk_collectmetrics {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/metrics", mode: "symlink"

    input:
        path ref_fasta
        path bwa_bqsr_bam
        path bt_bqsr_bam
    
    output:
        path "bwa/alignment_metrics.txt"
        path "bwa/insert_size_metrics.txt"
        path "bwa/insert_size_histogram.pdf"
        path "bowtie/alignment_metrics.txt"
        path "bowtie/insert_size_metrics.txt"
        path "bowtie/insert_size_histogram.pdf"

    script:
    """
    mkdir -p bwa bowtie

    echo "Collect metrics for BWA"

    gatk CollectAlignmentSummaryMetrics \
    -R ${ref_fasta} \
    -I ${bwa_bqsr_bam} \
    -O bwa/alignment_metrics.txt

    gatk CollectInsertSizeMetrics \
    INPUT=${bwa_bqsr_bam} \
    OUTPUT=bwa/insert_size_metrics.txt \
    HISTOGRAM_FILE=bwa/insert_size_histogram.pdf

    echo "Collect metrics for Bowtie"

    gatk CollectAlignmentSummaryMetrics \
    -R ${ref_fasta} \
    -I ${bt_bqsr_bam} \
    -O bowtie/alignment_metrics.txt

    gatk CollectInsertSizeMetrics \
    INPUT=${bt_bqsr_bam} \
    OUTPUT=bowtie/insert_size_metrics.txt \
    HISTOGRAM_FILE=bowtie/insert_size_histogram.pdf
    """

 }

 process aln_metrics {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/metrics", mode: "copy"

    input:
        path ref_fasta
        path markeddup_bam
        val aligner_name
    
    output:
        path "${aligner_name}/${markeddup_bam.baseName}_stats.txt"
        path "qualimap_bamqc/*"
        path "${aligner_name}/${aligner_name}_picard_metrics.txt"

    script:
    """
    echo "Running samtools stats"

    samtools stats ${markeddup_bam} > ${aligner_name}/${markeddup_bam.baseName}_stats.txt

    echo "Running qualimap"

    qualimap bamqc \
    -bam ${markeddup_bam} \
    -outdir ${aligner_name}/qualimap_bamqc \
    -outformat HTML \
    -nt 4

    echo "Runnning picard alignment summary metrics"

    picard CollectAlignmentSummaryMetrics \
    R=${ref_fasta} \
    I=${markeddup_bam} \
    O=${aligner_name}/${aligner_name}_picard_metrics.txt
    """
}

process gatk_haplotypecaller {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/gatk_hapcaller", mode: "symlink"

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path bwa_bqsr_bam
        val bwa_output_filename
    
    output:
        path "${bwa_output_filename}_raw_variants.vcf", emit: "bwa_raw_var"
    
    script:
    """
    echo "Run HaplotyeCaller on BWA"

    gatk HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${bwa_bqsr_bam} \
    -O ${bwa_output_filename}_raw_variants.vcf
    """
 }

 process gatk_selectsnpsindels {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/gatk_hapcaller", mode: "symlink"

    input:
        path ref_fasta
        path raw_variants
        val var_type1
        val var_type2
        val output_filename1
        val output_filename2
    
    output:
        path "${output_filename1}.vcf"
        path "${output_filename2}.vcf"

    script:
    """
    echo "Select variant types from BWA aligned sequence"

    gatk SelectVariants \
    -R ${ref_fasta} \
    -V ${raw_variants} \
    --select-type ${var_type1} \
    -O ${output_filename1}.vcf

    gatk SelectVariants \
    -R ${ref_fasta} \
    -V ${raw_variants} \
    --select-type ${var_type2} \
    -O ${output_filename2}.vcf
    """
 }

workflow {

    // Create input channel
    fasta_ch = Channel.fromPath(params.ref_fasta)
    read1_ch = Channel.fromPath(params.read1)
    read2_ch = Channel.fromPath(params.read2)
    refvcf_ch = Channel.fromPath(params.refvcf)

    // Create index file for input FASTA file
    st_index_out = samtools_index(fasta_ch)

    // Create sequence dictionary for input FASTA file
    gatk_dict_out = gatk_createseqdict(fasta_ch)

    // Create index file for reference VCF file
    gatk_vcfidx = gatk_indexvcf(refvcf_ch)

    // Run FASTQC to check read quality
    fastqc(read1_ch, read2_ch)

    // Remove adapter sequences from reads
    cutadapt_out = cutadapt(read1_ch, read2_ch)

    // Trim reads
    trimmed_out = trimmomatic(cutadapt_out.adapt_read1, cutadapt_out.adapt_read2)

    // Run BWA alignment
    bwa_out = BWA(fasta_ch, trimmed_out.trimmed_read1, trimmed_out.trimmed_read2, "SRR3340911")

    // Run Bowtie2 alignment
    bowtie_out = bowtie(fasta_ch, "tair10_ref", trimmed_out.trimmed_read1, trimmed_out.trimmed_read2, "SRR3340911")

    // Add missing RG lines to Bowtie .sam file
    rg_sam = add_rg(bowtie_out.bt_sam, "SRR3340911", "lib1", "ILLUMINA", "seq", "SRR3340911", "SRR3340911")

    // Convert SAM to BAM, sort and index 
    // st_sort_index()

    // Mark duplicates and sort the BWA sam file
    // markdupsort_out = markdupsort(bwa_out.bwa_sam, "bwa_SRR3340911", rg_sam.bt_rg_sam, "bowtie_SRR3340911")

    // Build BQSR model and apply it to the aligned reads file
    // bqsr_out = bqsr(fasta_ch, st_index_out.ref_fasta_fai, gatk_dict_out.ref_fasta_dict, refvcf_ch, gatk_vcfidx.ref_vcf_idx, markdupsort_out.bwa_sort_dedup_bam, "bwa_recal_data", markdupsort_out.bt_sort_dedup_bam, "bowtie_recal_data")

    // Collect alignment and insert size metrics
    // gatk_collectmetrics(fasta_ch, bqsr_out.bwa_bqsr_bam, bqsr_out.bt_bqsr_bam)

    // Alignment metrics
    // aln_metrics(fasta_ch, )

    // Variant calling via GATK HaplotypeCaller 
    // hapcaller_out = gatk_haplotypecaller(fasta_ch, st_index_out.ref_fasta_fai, gatk_dict_out.ref_fasta_dict, bqsr_out.bwa_bqsr_bam, "SRR3340911")

    // Selecting variant types
    // gatk_selectsnpsindels(fasta_ch, hapcaller_out.bwa_raw_var, "SNP", "INDEL", "raw_snps", "raw_indels")
}
