#!/usr/bin/env nextflow

// Predict m6A sites 

process FT_M6A {
    container 'docker://tamdbho/fiberseq:latest'
    label 'big_mem'
    tag "${id}"
    input:
        tuple val (id), path (bam_file)
    output:
        tuple val (id), path ("*_m6a.bam")
    script:
    def name = bam_file.getBaseName()
    """
    ft predict-m6a ${bam_file} ${name}_m6a.bam \
        -t 20 \
        -v \
        -k      
    """
}

// ALIGNMENT
process PBI_INDEX { 
    container 'docker://staphb/pbtk:3.1.1'
    publishDir "${params.out_dir}/${id}/m6a", mode: 'copy'
    input:
        tuple val (id), path (m6a_bam)
    output:
        tuple val (id), path(m6a_bam), path("*.pbi")
    script:
    """
    pbindex -j 10 ${m6a_bam}
    """
}

process ZMW_FILTER {
    container 'docker://staphb/pbtk:3.1.1'
    publishDir "${params.temp_dir}/${id}", mode: 'symlink'
    input:
        tuple val (id), path(m6a_bam), path(pbi)
    output:
        path "*_zmws.txt"
    script:
    def name = m6a_bam.getBaseName().replace('_m6a', '_zmws')
    """
    zmwfilter -j 10 --show-all ${m6a_bam} > "${name}.txt"
    """
}

process ZMW_INCLUDE {
    container 'docker://staphb/pbtk:3.1.1'
    publishDir "${params.temp_dir}/${id}", mode: 'symlink'
    input:
        tuple val (id), path(m6a_bam), path(pbi)
        path zmws_txt
    output:
        tuple val (id), path ("*_zmws.bam")
    script:
    def name = m6a_bam.getBaseName().replace('_m6a', '')
    """
    zmwfilter -j 10 --include ${zmws_txt} ${m6a_bam} ${name}_zmws.bam
    """  
}

process PBMM2_ALIGN {
    container 'docker://staphb/pbmm2:1.13.1'
    label 'big_mem'
    publishDir "${params.out_dir}/${id}/pbmm2", mode: 'copy'
    tag "${id}"
    input: 
        tuple val (id), path (zws_bam)
        path ref_file
	path ref_index
    output:
        tuple val (id),path ("*.bam") , path ("*.bai")
    
    script:
    def name = zws_bam.getBaseName().replace('_zmws', '_pbmm2')
    """
    pbmm2 align \
        -j 10 \
        --preset CCS --sort \
        --sort-memory 1G \
        --log-level INFO \
        --strip \
        --unmapped \
        ${ref_file} ${zws_bam} ${name}.bam
    """
}

//  BED FILES
process FT_EXTRACT {
    container 'docker://tamdbho/fiberseq:latest'
    publishDir "${params.out_dir}/${id}/ft_extract", mode: 'symlink'
    input:
        tuple val (id), path (pbmm2_bam), path (pbmm2_bai)
    output:
        tuple val (id), path ("*bed")
    script:
    def name = pbmm2_bam.getBaseName().replace('_pbmm2', '_ext')
    """
        ft extract ${pbmm2_bam} --m6a "${name}_m6a.bed" 
        ft extract ${pbmm2_bam} --cpg "${name}_cpg.bed" 
        ft extract ${pbmm2_bam} --msp "${name}_msp.bed" 
        ft extract ${pbmm2_bam} --nuc "${name}_nuc.bed" 
    """
}


process BED_SORT {
    container 'docker://staphb/bedtools:2.31.1'
    input:
        tuple val (id), path (bed_files)
    output:
        tuple val (id), path ("*_sorted.bed")
    script:
    def name = bed_files.getBaseName()
    """
    bedtools sort -i "${bed_files}" > "${name}_sorted.bed"
    """
}

process BED_INDEX {
    container 'docker://polumechanos/igvtools:latest'
    publishDir "${params.out_dir}/${id}/BED", mode: 'symlink'
    input:
        tuple val (id), path (sorted_bed)
    output:
        tuple val (id), path (sorted_bed), path ("*_sorted.bed.idx")
    script:
    """
    igvtools index ${sorted_bed}
    """
}

// SORT BAM FILES

process BAM_SORT {
    container 'docker://staphb/samtools:1.21'
    input:
        tuple val (id), path (pbmm2_bam), path (pbmm2_bai)
    output:
        tuple val (id), path ("*_sorted.bam")
    script:
    def name = pbmm2_bam.getBaseName().replace('_pbmm2', '_sorted')
    """
    samtools sort -@ 20 ${pbmm2_bam} -o ${name}.bam
    """
    
}

// MERGE BAM FILES

process BAM_MERGE {
    container 'docker://staphb/samtools:1.21'
    input:
        tuple val (id), path (sorted_bamList)
    output:
        tuple val (id), path ("*_merged.bam")
    script:
    def sorted_bamPair = sorted_bamList.collect { each -> each.toString() }.join(' ')
    """
    samtools merge -@ 20 -f ${sorted_bamPair} -o ${id}_merged.bam
    """
}

process BAM_INDEX {
    container 'docker://staphb/samtools:1.21'
    publishDir "${params.out_dir}/${id}/BAM", mode: 'copy'
    tag "${id}"
    input:
        tuple val (id), path (sorted_bam)
    output:
        tuple val (id), path (sorted_bam), path ("*.bai")
    script:
    """
    samtools index ${sorted_bam}
    """
}

// DEEPVARIANT

process DEEPVARIANT {
    container 'docker://google/deepvariant:1.8.0'
    publishDir "${params.temp_dir}/${id}", mode: 'symlink'
    label 'big_mem'
    input:
        tuple val (id), path (sorted_bam), path (sorted_bai)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*.vcf.gz"), path ("*.gvcf.gz")
    script:
    def name = sorted_bam.getBaseName().replace('_*', '_deepvariant')
    """
    run_deepvariant \
        --model_type="PACBIO" \
        --ref=${ref_file} \
        --reads=${sorted_bam} \
        --output_vcf=${id}_deepvariant.vcf.gz \
        --output_gvcf=${id}_deepvariant.gvcf.gz \
        --num_shards=8 \
        --sample_name=${id}
    """
}

process DEEPVARIANT_INDEX {
    container 'docker://quay.io/biocontainers/bcftools:1.7--0'
    publishDir "${params.out_dir}/${id}/deepvariant", mode: 'copy'
    input:
        tuple val (id), path (vcf), path (gvcf)
    output:
        tuple val (id), path (vcf), path ("*_deepvariant.vcf.gz.tbi"), emit: vcf
        tuple val (id), path (gvcf), path ("*_deepvariant.gvcf.gz.tbi"), emit: gvcf
    script:
    """
    bcftools index -t ${vcf}
    bcftools index -t ${gvcf}
    """
}

// PBSV

process PBSV_DISCOVER {
    container 'docker://quay.io/pacbio/pbsv:2.9.0_1.14_build1'
    input:
        tuple val (id), path (sorted_bam), path (sorted_bai)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*.svsig.gz")
    script:
    """
    pbsv discover \
    --ccs \
    --sample ${id} \
    ${sorted_bam} ${id}_pbsv.svsig.gz
    """
}

process PBSV {
    container 'docker://quay.io/pacbio/pbsv:2.9.0_1.14_build1'
    label 'big_mem'
    publishDir "${params.temp_dir}/${id}", mode: 'symlink'
    input:
        tuple val (id), path (sorted_bam), path (sorted_bai), path (svsig_file)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*pbsv.vcf")
    script:
    """
    pbsv call ${ref_file} \
    --ccs ${svsig_file} ${id}_pbsv.vcf
    """
}

process PBSV_INDEX {
    container 'docker://quay.io/pacbio/pbsv:2.9.0_1.14_build1'
    publishDir "${params.out_dir}/${id}/pbsv", mode: 'copy'
    input:
        tuple val (id), path (vcf)
    output:
        tuple val (id), path ("*pbsv.vcf.gz"), path ("*pbsv.vcf.gz.tbi")
    script:
    """
    bgzip -@ 4 -c ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

process SNIFFLES {
    container 'docker://quay.io/biocontainers/sniffles:2.5.3--pyhdfd78af_0'
    publishDir "${params.out_dir}/${id}/sniffles", mode: 'copy'
    label 'big_mem'
    input:
        tuple val (id), path (sorted_bam), path (sorted_bai)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*sniffles.vcf")
    script:
    """
    sniffles \
    --input ${sorted_bam} \
    --sample-id ${id} \
    --ref ${ref_file} \
    --vcf ${id}_sniffles.vcf
    """
}

process HIPHASE {
    container 'docker://quay.io/biocontainers/hiphase:1.4.5--h9ee0642_0'
    publishDir "${params.out_dir}/${id}/hiphase", mode: 'copy'
    label 'big_mem'
    input:
        tuple val (id), path (sorted_bam), path (sorted_bai), path(deepvariant_vcf), path(deepvariant_tbi), path(pbsv_vcf), path(pbsv_tbi)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*deepvariant.phased.vcf.gz") , path ("*deepvariant.phased.vcf.gz.tbi"), emit: deepvariant_vcf
        tuple val (id), path ("*pbsv.phased.vcf.gz") , path ("*pbsv.phased.vcf.gz.tbi"), emit: pbsv_vcf
        tuple val (id), path ("*_haplotagged.bam"), emit: bam
    script:
    """
    hiphase -t 16 \
    --ignore-read-groups \
    --bam ${sorted_bam} \
    --output-bam ${id}_haplotagged.bam \
    --reference ${ref_file} \
    --vcf ${deepvariant_vcf} \
    --output-vcf ${id}_deepvariant.phased.vcf.gz \
    --vcf ${pbsv_vcf} \
    --output-vcf ${id}_pbsv.phased.vcf.gz \
    --haplotag-file ${id}_read-level-phasing.tsv \
    --summary-file ${id}_summary.tsv \
    --stats-file ${id}_stats.tsv \
    --blocks-file ${id}_blocks.tsv 
    """
}

process FIRE {
    container 'docker://tamdbho/fiberseq:latest'
    publishDir "${params.out_dir}/${id}/FIRE", mode: 'copy'
    label 'big_mem'
    input:
        tuple val (id), path (haplotagged_bam)
        path ref_file
    output:
        tuple val (id), path ("*filtered.cram"), path ("*filtered.cram.crai")
    script:
    def threads = 32
    """
    samtools view -@ ${threads} -u -F ${params.filter_flag} ${haplotagged_bam} \
        | ft fire -F ${params.filter_flag} -t ${threads} \
            --min-msp ${params.min_msp} \
            --min-ave-msp-size ${params.min_ave_msp_size} \
            --skip-no-m6a \
            - - \
        | samtools view -C -@ ${threads} -T ${ref_file} \
            --output-fmt-option embed_ref=1 \
        | samtools view -C -@ ${threads} -T ${ref_file} \
            --output-fmt-option embed_ref=1 \
            --input-fmt-option required_fields=0x1bff \
            --write-index -o ${id}-fire-filtered.cram 
    """
}

// pb_CpG_calling
// using v3.0.0+: --model arugment is not needed
process pb_CpG_calling {
    container 'docker://tamdbho/pacbio:latest'
    publishDir "${params.out_dir}/${id}/CpG", mode: 'copy'
    label 'big_mem'
    input:
        tuple val(id), path (haplotagged_bam), path (haplotagged_bai)
        path (ref_file)
        path (ref_index)
    output:
        tuple val (id), path ("*.bed.gz"), emit: bed
        tuple val (id), path ("*.bw"), emit: bw
    script:
    def threads = 8
    """
    aligned_bam_to_cpg_scores \
        --bam ${haplotagged_bam} \
        --modsites-mode reference \
        --ref ${ref_file} \
        --output-prefix ${id}.cpg \
        --threads ${threads} \
        --pileup-mode model
    """
}