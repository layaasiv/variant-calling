/*
 * Pipeline parameters
 */

// Primary input
params.ref_fasta = "/home/layaasiv/professional-dev/variant-calling/nextflow/ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
params.read1 = "/home/layaasiv/professional-dev/variant-calling/reads/SRR3340911_1.fastq"
params.read2 = "/home/layaasiv/professional-dev/variant-calling/reads/SRR3340911_2.fastq"
params.refvcf = "/home/layaasiv/professional-dev/variant-calling/reference/arabidopsis_thaliana.vcf"
params.truth_vcf = "/home/layaasiv/professional-dev/variant-calling/var_caller/truth_renamed.vcf.gz"

// Output directory
params.outdir = "/home/layaasiv/professional-dev/variant-calling/nf_revised/results_nf"

/*
 * Generate fasta index file
 */

process samtools_index {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

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

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

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

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

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

    publishDir "${params.outdir}/preprocessing", mode: "symlink"

    input:
        path input_read1
        path input_read2

    output:
        path "fastqc_reports/*"
    
    script:
    """
    # Create the output directory if it doesn't exist
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
        path "adapt_${input_read1.baseName}", emit: "adapt_read1"
        path "adapt_${input_read2.baseName}", emit: "adapt_read2"
    
    script:
    """
    cutadapt \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o adapt_${input_read1.baseName} \
    -p adapt_${input_read2.baseName} \
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
        tuple path("*.sam"), val(base_fn), emit: "bwa_sam"
    
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

    publishDir "${params.outdir}/bowtie2", mode: "symlink"

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
        path "${base_fn}_bowtie2_raw.sam", emit: "raw_bt_sam"

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
    ${base_fn}_bowtie2_raw.sam
    """
}

/*
 Adding read group information to the Bowtie2 .sam file
 */

process add_rg {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/bowtie2", mode: "symlink"

    input:
        path bowtie_aligned_sam
        val RGID
        val RGLB
        val RGPL
        val RGPU
        val RGSM
        val base_fn
    
    output:
        tuple path("*.sam"), val(base_fn), emit: "bowtie2_sam"

    script:
    """
    picard AddOrReplaceReadGroups \
    I=${bowtie_aligned_sam} \
    O=${base_fn}_bowtie2.sam \
    RGID=${RGID} \
    RGLB=${RGLB} \
    RGPL=${RGPL} \
    RGPU=${RGPU} \
    RGSM=${RGSM}
    """
}

/*
Convert SAM to BAM, sort, and index files
*/

process st_sort_index {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    input:
        tuple path(raw_aligned_sam), val(base_fn)
        val aligner

    publishDir "${params.outdir}/${aligner}", mode: "symlink"

    output:
        tuple path("${base_fn}_${aligner}_sorted.bam"), path("${base_fn}_${aligner}_sorted.bam.bai"), val(base_fn), val(aligner), emit: "sorted_bam_bai"

    script:
    """
    samtools view -Sb ${raw_aligned_sam} | samtools sort -o ${base_fn}_${aligner}_sorted.bam

    samtools index ${base_fn}_${aligner}_sorted.bam
    """
}

/*
Mark PCR duplicates using picard 
*/

process markduplicates {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/${aligner}", mode: "symlink"

    input:
        tuple path(sorted_bam), path(sorted_bai), val(base_fn), val(aligner)

    output:
        tuple path("${base_fn}_${aligner}_markeddup.bam"), path("${base_fn}_${aligner}_markeddup.bai"), val(base_fn), val(aligner), emit: "markdup_bam_bai"
        path "${base_fn}_${aligner}_markdup_metrics.txt"

    script:
    """
    picard MarkDuplicates \
    I=${sorted_bam} \
    O=${base_fn}_${aligner}_markeddup.bam \
    M=${base_fn}_${aligner}_markdup_metrics.txt \
    CREATE_INDEX=true
    """
}

/*
Collect alignment metrics
*/

process alignment_metrics {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    input:
        path ref_fasta
        tuple path(markeddup_bam), path(markeddup_bai), val (base_fn), val(aligner)

    publishDir "${params.outdir}/metrics/aligners/${aligner}", mode: "copy"
    
    output:
        path "${base_fn}_${aligner}_st_stats.txt"
        path "qualimap_bamqc/*"
        path "${base_fn}_${aligner}_picard_metrics.txt"
    
    script:
    """
    echo "Running samtools stats"

    samtools stats ${markeddup_bam} > ${base_fn}_${aligner}_st_stats.txt

    echo "Running qualimap"

    qualimap bamqc \
    -bam ${markeddup_bam} \
    -outdir qualimap_bamqc \
    -outformat HTML \
    -nt 4

    echo "Runnning picard alignment summary metrics"

    picard CollectAlignmentSummaryMetrics \
    R=${ref_fasta} \
    I=${markeddup_bam} \
    O=${base_fn}_${aligner}_picard_metrics.txt
    """
}

/*
 * Base quality score recalibration
 */

process bqsr {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/bwa", mode: "symlink"

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path ref_vcf
        path ref_vcf_idx
        path bwa_sorted_dedup_bam
        path bwa_sorted_dedup_bai
        val base_fn

    output:
        path "*.table", emit: "bwa_bqsr_recal_table"
        path "*_bqsr.bam", emit: "bwa_bqsr_bam"

    script:
    """
    echo "Building the BQSR model"

    gatk BaseRecalibrator \
    -I ${bwa_sorted_dedup_bam} \
    -R ${ref_fasta} \
    --known-sites ${ref_vcf} \
    -O ${base_fn}_bqsr.table

    echo "Applying the BQSR model"

    gatk ApplyBQSR \
    -I ${bwa_sorted_dedup_bam} \
    -R ${ref_fasta} \
    --bqsr-recal-file ${base_fn}_bqsr.table \
    -O ${base_fn}_bqsr.bam
    """
}

/*
* GATK HaplotypeCaller for variant calling
*/

process gatk_haplotypecaller {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/gatk_hapcaller", mode: "symlink"

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path bwa_bqsr_bam
        val base_fn
    
    output:
        path "${base_fn}_raw_var_HC.vcf", emit: "bwa_hapcaller_vcf"
    
    script:
    """
    echo "Run HaplotyeCaller"

    gatk HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${bwa_bqsr_bam} \
    -O ${base_fn}_raw_variants.vcf
    """
 }

/*
* FreeBayes variant caller
*/

process freebayes {
    conda "/home/layaasiv/miniconda3/envs/sra-tools"

    publishDir "${params.outdir}/freebayes", mode: "symlink"

    input:
        path ref_fasta
        path bwa_bqsr_bam
        val base_fn

    output:
        path "${base_fn}_raw_var_FB.vcf", emit: "bwa_fb_vcf"
    
    script:
    """
    freebayes -f ${ref_fasta} \
    ${bwa_bqsr_bam} > \
    ${base_fn}_raw_var_FB.vcf

    """
}

/*
* Evaluating variant callers: GATK VariantEval
*/

process gatk_varianteval {
    container "broadinstitute/gatk"

    publishDir "${params.outdir}/metrics/var-callers", mode: "copy"

    input:
        path ref_fasta
        path ref_fai
        path ref_fasta_dict
        path query_vcf
        path query_vcf_idx
        path truth_vcf
        path truth_vcf_idx
        val var_caller
        val base_fn

    output:
        path "${base_fn}_${var_caller}_gatkvareval_output.eval"

    script:
    """
    gatk VariantEval \
    -R ${ref_fasta} \
    -eval ${query_vcf} \
    -comp ${truth_vcf} \
    -O ${base_fn}_${var_caller}_gatkvareval_output.eval
    """
}

/*
Evaluating variant callers: bcftools isec
*/ 

process bcftools_eval_gen {
    conda "/home/layaasiv/miniconda3/envs/happy-tools"

    publishDir "${params.outdir}/metrics/var-callers", mode: "copy"

    input:
        path query_vcf1
        path query_vcf1_idx
        path query_vcf2
        path query_vcf2_idx
        path truth_vcf
        val base_fn

    output:
        path "${base_fn}_truth.stats"
        path "isec_output/*"

    script:
    """
    bcftools isec -p isec_output -n=2 -w1 ${query_vcf1} ${query_vcf2}

    bcftools stats ${truth_vcf} > ${base_fn}_truth.stats
    """
}

/*
Evaluating variant callers: bcftools stats
*/ 

process bcftools_eval_vc_spec {
    conda "/home/layaasiv/miniconda3/envs/happy-tools"

    publishDir "${params.outdir}/metrics/var-callers", mode: "copy"

    input:
        path query_vcf
        path query_vcf_idx
        val base_fn
        val var_caller

    output:
        path "${base_fn}_${var_caller}.stats"

    script:
    """
    bcftools stats ${query_vcf} > ${base_fn}_${var_caller}.stats
    """
}

/*
Evaluating variant callers: Illumina hap.py
*/ 

process hap_py_eval {
    container "miguelpmachado/hap.py_illumina:0.3.10-8-gdac84d9-02"

    publishDir "${params.outdir}/metrics/var-callers", mode: "copy"

    input:
        path ref_fasta
        path ref_fai
        path ref_fasta_dict
        path truth_vcf
        path truth_vcf_idx
        path query_vcf
        path query_vcf_idx
        val var_caller

    output:
        path "happy/happy_${var_caller}*"

    script:
    """
    mkdir -p happy

    hap.py \
    ${truth_vcf} ${query_vcf} \
    -r ${ref_fasta} \
    -o happy/happy_${var_caller}
    """
}

workflow preprocess_qc {
    take:
    ref_fasta
    read1
    read2
    ref_vcf

    main:
    // Create index file for input FASTA file
    st_index_out = samtools_index(ref_fasta)
    // Create sequence dictionary for input FASTA file
    gatk_dict_out = gatk_createseqdict(ref_fasta)
    // Create index file for reference VCF file
    gatk_vcfidx = gatk_indexvcf(ref_vcf)
    // Run FASTQC to check read quality
    fastqc(read1, read2)
    // Remove adapter sequences from reads
    cutadapt_out = cutadapt(read1, read2)
    // Trim reads
    trimmed_out = trimmomatic(cutadapt_out.adapt_read1, cutadapt_out.adapt_read2)

    emit:
    trimmed_read1 = trimmed_out.trimmed_read1
    trimmed_read2 = trimmed_out.trimmed_read2
    ref_fasta_fai = st_index_out.ref_fasta_fai
    ref_fasta_dict = gatk_dict_out.ref_fasta_dict
    ref_vcf_idx = gatk_vcfidx.ref_vcf_idx

}

workflow bwa_wf {
    take:
    ref_fasta
    trimmed_read1
    trimmed_read2
    base_fn
    aligner

    main:
    // Run BWA aligner
    bwa_out = BWA(ref_fasta, trimmed_read1, trimmed_read2, base_fn)

    // Convert SAM to BAM, sort, and index
    st_sort_index_out = st_sort_index(bwa_out.bwa_sam, aligner)
    // Mark PCR duplicates
    markdup_out = markduplicates(st_sort_index_out.sorted_bam_bai)
    // Collect alignment metrics for BWA alignment
    alignment_metrics(ref_fasta, markdup_out.markdup_bam_bai)

    emit:
    bwa_markdup_bam = markdup_out.markdup_bam_bai.map { it[0] }
    bwa_markdup_bai = markdup_out.markdup_bam_bai.map { it[1] }

}


workflow bowtie2_wf {
    take:
    ref_fasta
    trimmed_read1
    trimmed_read2
    base_fn
    aligner
    idx_filename_base
    RGID
    RGLB
    RGPL
    RGPU
    RGSM

    main:
    // Run bowtie2 aligner
    bowtie2_out = bowtie(ref_fasta, idx_filename_base, trimmed_read1, trimmed_read2, base_fn)
    // Add read groups to the bowtie2 SAM
    bowtie2_rg_out = add_rg(bowtie2_out.raw_bt_sam, RGID, RGLB, RGPL, RGPU, RGSM, base_fn)

    // Convert to BAM, sort, and index
    st_sort_index_out = st_sort_index(bowtie2_rg_out.bowtie2_sam, aligner)
    // Mark PCR duplicates
    markdup_out = markduplicates(st_sort_index_out.sorted_bam_bai)
    // Collect alignment metrics for bowtie2 alignment
    alignment_metrics(ref_fasta, markdup_out.markdup_bam_bai)
}

workflow var_call_wf {
    take:
    ref_fasta
    ref_fasta_fai
    ref_fasta_dict
    ref_vcf
    ref_vcf_idx
    bwa_markdup_bam
    bwa_markdup_bai
    base_fn

    main:
    // Run GATK base quality recalibrator 
    bqsr_out = bqsr(ref_fasta, ref_fasta_fai, ref_fasta_dict, ref_vcf, ref_vcf_idx, bwa_markdup_bam, bwa_markdup_bai, base_fn)

    // Run GATK HaplotypeCaller
    hapcaller_out = gatk_haplotypecaller(ref_fasta, ref_fasta_fai, ref_fasta_dict, bqsr_out.bwa_bqsr_bam, base_fn)

    // Run FreeBayes variant caller
    freebayes_out = freebayes(ref_fasta, bqsr_out.bwa_bqsr_bam, base_fn)

    emit:
    hc_vcf = hapcaller_out.bwa_hapcaller_vcf
    fb_vcf = freebayes_out.bwa_fb_vcf

} 

workflow vc_eval_general_setup {
    take:
    truth_vcf

    main:
    // Index feature file for truth
    gatk_indexvcf_out = gatk_indexvcf(truth_vcf)

    emit:
    truth_vcf_idx = gatk_indexvcf_out.vcf_idx
    
}

workflow hc_eval {
    take:
    ref_fasta
    ref_fai
    ref_fasta_dict
    hc_vcf
    truth_vcf
    truth_vcf_idx
    hc_var_caller
    base_fn

    main:
    // Index feature file for query
    gatk_indexvcf_out = gatk_indexvcf(hc_vcf)
    
    // VC evaluation: GATK VariantEval
    gatk_varianteval_out = gatk_varianteval(ref_fasta, ref_fai, ref_fasta_dict, hc_vcf, gatk_indexvcf_out.vcf_idx, truth_vcf, truth_vcf_idx, hc_var_caller, base_fn)

    //VC evaluation: BCFTools stats 
    bcftools_eval_vc_spec(hc_vcf, gatk_indexvcf_out.vcf_idx, base_fn, hc_var_caller)

    // VC evaluation: hap.py
    hap_py_eval(ref_fasta, ref_fai, ref_fasta_dict, truth_vcf, truth_vcf_idx, hc_vcf, gatk_indexvcf_out.vcf_idx, hc_var_caller)

    emit:
    hc_vcf_idx = gatk_indexvcf_out.vcf_idx

}

workflow fb_eval {
    take:
    ref_fasta
    ref_fai
    ref_fasta_dict
    fb_vcf
    truth_vcf
    truth_vcf_idx
    fb_var_caller
    base_fn

    main:
    // Index feature file for query
    gatk_indexvcf_out = gatk_indexvcf(fb_vcf)
    
    // VC evaluation: GATK VariantEval
    gatk_varianteval_out = gatk_varianteval(ref_fasta, ref_fai, ref_fasta_dict, fb_vcf, gatk_indexvcf_out.vcf_idx, truth_vcf, truth_vcf_idx, fb_var_caller, base_fn)

    //VC evaluation: BCFTools stats 
    bcftools_eval_vc_spec(fb_vcf, gatk_indexvcf_out.vcf_idx, base_fn, fb_var_caller)

    // VC evaluation: hap.py
    hap_py_eval(ref_fasta, ref_fai, ref_fasta_dict, truth_vcf, truth_vcf_idx, fb_vcf, gatk_indexvcf_out.vcf_idx, fb_var_caller)

    emit:
    fb_vcf_idx = gatk_indexvcf_out.vcf_idx
}

workflow overall_eval {
    take:
    query_vcf1
    query_vcf1_idx
    query_vcf2
    query_vcf2_idx
    truth_vcf
    base_fn

    main:
    // VC evaluation: BCFTools isec & stats 
    bcftools_eval_gen(query_vcf1, query_vcf1_idx, query_vcf2, query_vcf2_idx, truth_vcf, base_fn)
}


workflow {
    // Raw inputs
    fasta_ch = Channel.fromPath(params.ref_fasta)
    read1_ch = Channel.fromPath(params.read1)
    read2_ch = Channel.fromPath(params.read2)
    refvcf_ch = Channel.fromPath(params.refvcf)

    // Alignment inputs
    base_fn_ch = Channel.value("SRR3340911")
    aligner_bwa_ch = Channel.value("bwa")
    aligner_bowtie_ch = Channel.value("bowtie2")

    // Bowtie2-specific inputs
    idx_base_fn_ch = Channel.value("tair10_ref")
    RGID_ch = Channel.value("SRR3340911")
    RGLB_ch = Channel.value("lib1")
    RGPL_ch = Channel.value("ILLUMINA")
    RGPU_ch = Channel.value("seq")
    RGSM_ch = Channel.value("SRR3340911")

    // Variant calling inputs
    truth_vcf_ch = Channel.fromPath(params.truth_vcf)
    hc_var_caller_ch = Channel.value("HC")
    fb_var_caller_ch = Channel.value("FB")

    // Preprocess data and complete QC of raw FASTQ read files
    preprocess_qc(
        fasta_ch, 
        read1_ch, 
        read2_ch, 
        refvcf_ch
        )

    bwa_wf(
        fasta_ch, 
        preprocess_qc.out.trimmed_read1, 
        preprocess_qc.out.trimmed_read2, 
        base_fn_ch, 
        aligner_bwa_ch
        )

    bowtie2_wf(
        fasta_ch, 
        preprocess_qc.out.trimmed_read1, 
        preprocess_qc.out.trimmed_read2, 
        base_fn_ch, 
        aligner_bowtie_ch, 
        idx_base_fn_ch,
        RGID_ch, 
        RGLB_ch,
        RGPL_ch,
        RGPU_ch,
        RGSM_ch
        )
    
    var_call_wf(
        fasta_ch,
        preprocess_qc.out.ref_fasta_fai,
        preprocess_qc.out.ref_fasta_dict,
        refvcf_ch,
        preprocess_qc.out.ref_vcf_idx,
        bwa_wf.out.bwa_markdup_bam,
        bwa_wf.out.bwa_markdup_bai,
        base_fn_ch
        )

    vc_eval_general_setup(
        truth_vcf_ch
        )

    hc_eval(
        fasta_ch,
        preprocess_qc.out.ref_fasta_fai,
        preprocess_qc.out.ref_fasta_dict,
        var_call_wf.out.hc_vcf,
        truth_vcf_ch,
        vc_eval_general_setup.out.truth_vcf_idx,
        hc_var_caller_ch,
        base_fn_ch
        )

    fb_eval(
        fasta_ch,
        preprocess_qc.out.ref_fasta_fai,
        preprocess_qc.out.ref_fasta_dict,
        var_call_wf.out.fb_vcf,
        truth_vcf_ch,
        vc_eval_general_setup.out.truth_vcf_idx,
        fb_var_caller_ch,
        base_fn_ch
        )

    overall_eval(
        var_call_wf.out.hc_vcf,
        hc_eval.out.hc_vcf_idx,
        var_call_wf.out.fb_vcf,
        fb_eval.out.fb_vcf_idx,
        truth_vcf_ch,
        base_fn_ch
        )
}