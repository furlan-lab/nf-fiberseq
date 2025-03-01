#!/usr/bin/env nextflow

process GET_CHROM {
    input:
        path ref_file
        path ref_index
    output:
        path 'chromosomes.txt'
    script:
    """
    get_chroms.py \
        --input ${ref_file} \
        --min_contig_length ${params.min_contig_length} > chromosomes.txt
    """
}

// FIRE OUTPUTS
process FIRE_SITES_CHROM {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire", mode: 'symlink'
    label 'big_mem'
    tag "${id}-${chromosome}"
    input:
        tuple val (id), path (cram), path (crai), val (chromosome)
    output:
        tuple val (id), path ("*-${chromosome}-sorted.bed.gz")
    script:
    def threads = 4
    """
    samtools view -@ ${threads} -u ${cram} ${chromosome} \
    | ft fire -t ${threads} --extract - \
        | LC_ALL=C sort --parallel=${threads} \
            -k1,1 -k2,2n -k3,3n -k4,4 \
        | bioawk -tc hdr '\$10<=${params.min_fire_fdr}' \
        | (grep '\\S' || true) \
        | (grep -v '^#' || true) \
        | bgzip -@ ${threads} \
    > ${id}-${chromosome}-sorted.bed.gz
    """
}


process FIRE_SITES_INDEX {
    container 'docker://quay.io/pacbio/pbsv:2.9.0_1.14_build1'
    publishDir "${params.out_dir}/${id}/FIRE/fire-sites", mode: 'copy'
    input:
        tuple val (id), path (fire_sites)
    output:
        tuple val (id), path ("*-fire-elements.bed.gz"), path ("*fire-elements.bed.gz.tbi")
    script:
    def fire_beds = fire_sites.collect { each -> each.toString() }.join(' ')
    """
    cat ${fire_beds} > "${id}-fire-elements.bed.gz"
    tabix -p bed ${id}-fire-elements.bed.gz
    
    """
}

// COVERAGE INFORMATION

process BEDGRAPH {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire", mode: 'symlink'
    label 'big_mem'
    input:
        tuple val (id), path (cram), path (crai)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*bed.gz"), path ("*bed.gz.tbi")
    script:
    """
    mosdepth -f ${ref_file} -t 16 tmp ${cram}
    bgzip -cd tmp.per-base.bed.gz \
            | LC_ALL=C sort --parallel=16 -k1,1 -k2,2n -k3,3n -k4,4  \
            | bgzip -@ 16 \
        > ${id}.bed.gz
    tabix -f -p bed ${id}.bed.gz
    """
}

process COVERAGE {
    publishDir "${params.out_dir}/${id}/FIRE/coverage", mode: 'copy' 
    input:
        tuple val(id), path(bedgraph), path(bedgraph_tbi)
        path chromosomes
    output:
        tuple val (id), path ("*median-coverage.txt") , emit: median
        tuple val (id), path ("*minimum-coverage.txt") , emit: minimum
        tuple val (id), path ("*maximum-coverage.txt") , emit: maximum
    script:
    """
    cov.py \
        --input ${bedgraph} \
        --output_median ${id}-median-coverage.txt \
        --output_min ${id}-minimum-coverage.txt \
        --output_max ${id}-maximum-coverage.txt \\
        --coverage_within_n_sd ${params.coverage_within_n_sd} \
        --min_coverage ${params.min_coverage} \
        --chromosomes ${chromosomes}
    """
}

process FIBER_LOCATIONS_CHROM {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire", mode: 'symlink'
    label 'big_mem'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), path(cram), path(crai), val (chromosome)
    output:
        tuple val (id), path ("*${chrom}.fiber-locations.bed.gz")
    script:
    def threads = 4
    """
    (samtools view -@ ${threads} -u ${cram} ${chromosome} \
        | ft extract -t ${threads} -s --all - \
        | hck -F '#ct' -F st -F en -F fiber -F strand -F HP ) \
        | (grep -v "^#" || true) \
        | bgzip -@ ${threads} \
    > ${id}-${chrom}.fiber-locations.bed.gz
    """

}

process FIBER_LOCATIONS_INDEX {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire", mode: 'symlink'
    label 'big_mem'
    input:
        tuple val(id), path(fibers_beds), path (bedgraph), path (bedgraph_tbi), path (min_cov), path (max_cov)
    output:
        tuple val (id), path ("*-filtered-fiber-locations.bed.gz"), path ("*-filtered-fiber-locations.bed.gz.tbi"), emit: filtered
        tuple val (id), path ("*-fiber-location.bed.gz"), path ("*-fiber-location.bed.gz.tbi"), emit: bed
    script:
    def threads = 4
    def beds = fibers_beds.collect { each -> each.toString() }.join(' ')
    """
    cat ${beds} > ${id}-fiber-location.bed.gz
    tabix -f -p bed ${id}-fiber-location.bed.gz

    MIN=$(cat ${min_cov})
    MAX=$(cat ${max_cov})
    bedtools intersect -header -sorted -v -f ${params.max_frac_overlap} \
        -a ${id}-fiber-locations.bed.gz \
        -b <(bgzip -cd ${bedgraph} | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX') \
    | bgzip -@ ${threads} \
    > ${id}-filtered-fiber-locations.bed.gz
    tabix -f -p bed ${id}-filtered-fiber-locations.bed.gz
    """
}
// FIRE PEAKS


// TRACKHUB
// this function filters out reads with 0 coverage and also reads from some annotations beds
process SHUFFLE_EXCLUDE {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.out_dir}/${id}/FIRE/exclude", mode: 'copy'
    input: 
        tuple val(id), path(fiber_locations), path(fiber_locations_tbi)
        path excludes
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*exclude-from-shuffles.bed.gz")
    script:
    def threads = 4
    """
    ( \
        bedtools genomecov -bga -i ${fiber_locations} -g ${ref_index} | awk '\$4 == 0'; \
        less ${excludes} \
    ) \
        | cut -f 1-3 \
        | bedtools sort \
        | bedtools merge \
        | bgzip -@ ${threads} \
    > ${id}-exclude-from-shuffles.bed.gz
    """

}

process FILTERED_SHUFFLED_FIBER_LOCATIONS_CHROM {
    // add bioawk and use pacbio docker for this one so we don't have to use big_mem
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire/shuffled", mode: 'symlink'
    label 'big_mem'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), path(fiber_locations), path(fiber_locations_tbi), path (exclude_bed), val (chromosome)
        path ref_file
        path ref_index
    output:
        tuple val (id), val (chromosome), path ("*${chromosome}.fiber-locations-shuffled.bed.gz")
    script:
    def threads = 4
    """
    tabix ${fiber_locations} ${chromosome} \
        | bioawk -t '{print $1,$2,$3,$4,$2}' \
        | bedtools shuffle -chrom -seed 42 \
            -excl ${exclude_bed} \
            -i - \
            -g ${ref_index} \
        |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
        | bgzip -@ ${threads} \
    > ${id}-${chromosome}.fiber-locations-shuffled.bed.gz
    """
}

process SHUFFLED_PILEUP_CHROM {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire/pileup", mode: 'symlink'
    label 'big_mem'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), path(cram), path(crai), val (chromosome), path (filtered_shuffled_beds) 
    output:
        tuple val (id), path ("*${chromosome}.-shuffled.pileup.bed.gz")
    script:
    def threads = 4
    """
    ft pileup ${cram} ${chromosome} -t ${threads} \
        --fiber-coverage --shuffle ${filtered_shuffled_beds} \
        --no-msp --no-nuc \
        | bgzip -@ ${threads} \
    > ${id}-${chromosome}-shuffled.pileup.bed.gz
    """
}

process SHUFFLED_PILEUP {
    container 'docker://tamdbho/fire:latest'
    publishDir "${params.out_dir}/${id}/FIRE/pileup", mode: 'copy'
    input:
        tuple val(id), path (pileup_chroms)
    output:
        tuple val (id), path ("*shuffled-pileup.bed.gz"), path ("*shuffled-pileup.bed.gz.tbi")
    script:
    def threads = 4
    def pileup_chrom = pileup_chroms.collect { each -> each.toString() }.join(' ')
    """
    cat ${pileup_chrom} > ${id}-shuffled-pileup.bed.gz
    tabix -f -p bed ${id}-shuffled-pileup.bed.gz
    """
}

process FDR_TABLE {
    publishDir "${params.out_dir}/${id}/FIRE/fire-peaks", mode: 'copy'
    input:
        tuple val(id), path (shuffle_pileup), path (shuffle_pileup_tbi), path (min_cov), path (max_cov)
    output:
        tuple val (id), path ("*-fire-score-to-fdr.tbl")
    script:
    """
    MIN=$(cat ${min_cov})
    MAX=$(cat ${max_cov})
    fdr-table.py \
    -v 1 ${shuffle_pileup} ${id}-*-fire-score-to-fdr.tbl --max-cov $MAX --min-cov $MIN
    """
}

process PILEUP_CHROM {
    container 'docker://tamdbho/fireseq:latest'
    publishDir "${params.temp_dir}/${id}/fire/pileup", mode: 'symlink'
    label 'big_mem'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), path (cram), path (crai), val (chromosome)
    output:
        tuple val (id), val (chromosome), path ("*${chromosome}.pileup.bed.gz")
    script:
    def threads = 4
    """
    ft pileup -t ${threads} \
        --haps --fiber-coverage \
        ${cram} ${chromosome} \
        | bgzip -@ ${threads} \
        > ${id}-${chromosome}.pileup.bed.gz
    """
}

process FDR_TRACK_CHROM {
    publishDir "${params.temp_dir}/${id}/fire/pileup", mode: 'symlink'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), val (chromosome), path (pileup_chroms), path (fdr_table)
    output:
        tuple val (id), path ("*${chromosome}.FDR.track.bed")
    script:
    """
    fdr-table.py -v 1 \
    --fdr-table ${fdr_tbl} \
    ${pileup_chroms} ${id}-${chromosome}.FDR.track.bed
    """
}

process PILEUP {
    publishDir "${params.temp_dir}/${id}/fire/pileup", mode: 'symlink'
    input:
        tuple val(id), path (fdr_pileup_beds)
    output:
        tuple val (id), path ("*pileup.bed.gz"), path ("*pileup.bed.gz.tbi") , emit: bed
        tuple val (id), path ("*pileup.fofn") , emit: fofn
    script:
    def beds = fdr_pileup_beds.collect { each -> each.toString() }.join(' ')
    def threads = 4
    """
    printf '\nMaking FOFN\n'
    echo ${beds} > ${id}-pileup.fofn  
    printf '\nMake header\n'
    ((cat $(cat ${id}-pileup.fofn) | grep "^#" | head -n 1) || true) \
        | bgzip -@ ${threads} \
        > ${id}-pileup.bed.gz

    printf '\nConcatenating\n'
    cat $(cat ${id}-pileup.fofn \
        | grep -v "^#" \
        | bgzip -@ ${threads} \
    >> ${id}-pileup.bed.gz

    printf '\nIndexing\n'
    tabix -f -p bed ${id}-pileup.bed.gz
    """
}

process HELPER_FDR_PEAKS_BY_FIRE_ELEMENTS {
    container 'docker://tamdbho/fire:latest'
    publishDir "${params.temp_dir}/${id}/fire/fire-peaks", mode: 'symlink'
    tag "${id}-${chromosome}"
    input:
        tuple val(id), path (pileup_bed), path (pileup_tbi), path (firesite_bed), path (firesite_tbi), val (chromosome)
    output:
        tuple val (id), val(chromosome), path ("*${chromosome}-fire-peaks.bed.gz")
    script:
    def threads = 2
    """
    HEADER=$(bgzip -cd ${pileup_bed} | head -n 1 || true)
    NC=$(echo $HEADER | awk '{{print NF}}' || true)
    FIRE_CT=$((NC+1))
    FIRE_ST=$((NC+2))
    FIRE_EN=$((NC+3))
    FIRE_SIZE=$((NC+4))
    FIRE_ID=$((NC+5))

    OUT_HEADER=$(printf "$HEADER\\tpeak_chrom\\tpeak_start\\tpeak_end\\tFIRE_IDs\\tFIRE_size_mean\\tFIRE_size_ssd\\tFIRE_start_ssd\\tFIRE_end_ssd")
    echo $OUT_HEADER
        
    ( \
        printf "$OUT_HEADER\\n"; \
        tabix -h ${pileup_bed} ${chromosome} \
            | bioawk -tc hdr '(NR==1)||($is_local_max=="true")' \
            | csvtk filter -tT -C '$' -f "FDR<=${params.max_peak_fdr}" \
            | csvtk filter -tT -C '$' -f "fire_coverage>1" \
            | bioawk -tc hdr '(NR==1)||($fire_coverage/$coverage>=${params.min_per_acc_peak})' \
            | bedtools intersect -wa -wb -sorted -a - \
                -b <(tabix ${firesite_bed} ${chromosome} \
                        | cut -f 1-3 \
                        | awk -v OFMT="%f" '{print $0"\t"$3-$2"\t"NR}' \
                    ) \
            | bedtools groupby -g 1-$NC \
                -o first,median,median,collapse,mean,sstdev,sstdev,sstdev \
                -c $FIRE_CT,$FIRE_ST,$FIRE_EN,$FIRE_ID,$FIRE_SIZE,$FIRE_SIZE,$FIRE_ST,$FIRE_EN \
    ) \
        | hck -f 1,$FIRE_ST,$FIRE_EN,2-$NC,$FIRE_SIZE- \
        | csvtk round -tT -C '$' -n 0 -f 2,3 \
        | bedtools sort -header -i - \
        | bgzip -@ ${threads} \
        > ${id}-${chromosome}-fire-peaks.bed.gz
    """
}

process FDR_PEAKS_BY_FIRE_ELEMENTS_CHROM {
    container 'docker://tamdbho/fire:latest'
    publishDir "${params.temp_dir}/${id}/fire/fire-peaks", mode: 'symlink'
    tag "${id}-${chromosome}"
    input:
        tuple val (id), val(chromosome), path (firepeak_bed), path (min_cov), path (max_cov)
    output:
        tuple val (id), path ("*${chromosome}-grouped-fire-peaks.bed.gz")
    script:
    def threads = 4
    """
    bgzip -cd ${firepeak_bed} \
        | merge_fire_peaks.py -v 1 \
            --max-cov $(cat {max_cov}) \
            --min-cov $(cat {min_cov}) \
            --min-frac-accessible ${params.min_frac_accessible} \
        | bgzip -@ ${threads} \
    > ${id}-${chromosome}-grouped-fire-peaks.bed.gz
    """
}



process FIRE_PEAKS {
    container 'docker://tamdbho/fire:latest'
    publishDir "${params.out_dir}/${id}/FIRE/fire-peaks", mode: 'copy'
    input:
        tuple val (id), path (grouped_firepeak_bed)
    output:
        tuple val (id), path ("*-fire-peaks.bed.gz"), path ("*fire-peaks.bed.gz.tbi"), emit: bed
        tuple val (id), path ("*-fire-peaks.fofn"), emit: fofn
    script:
    def beds = grouped_firepeak_bed.collect { each -> each.toString() }.join(' ')
    def threads = 4
    """
    printf "\nMaking FOFN\n"
    echo ${beds} > ${id}-fire-peaks.fofn

    printf "\nMaking header\n"        
    ((cat $(cat ${id}-fire-peaks.fofn) | bgzip -cd | grep "^#" | head -n 1) || true) \
        | bgzip -@ ${threads} > ${id}-fire-peaks.bed.gz

    printf "\nConcatenating\n"
    cat $(cat ${id}-fire-peaks.fofn)) | bgzip -cd -@ ${threads} | grep -v "^#" \
        | bgzip -@ ${threads} >> ${id}-fire-peaks.bed.gz
    
    printf "\nIndexing\n"
    tabix -f -p bed ${id}-fire-peaks.bed.gz
"""
}

// TRACKHUB
process FIRE_PEAKS_BB {
    container 'docker://tamdbho/fire:latest'
    publishDir "${params.out_dir}/${id}/trackHub/bb/fire-peaks.bb", mode: 'copy'
    input:
        tuple val (id), path (firepeak_bed)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*-fire-peaks.bb")
    script:
    def template = "./template.fire_peak.as"
    """
    bgzip -cd ${firepeak_bed} \
        | bioawk -tc hdr '{{print $1,$2,$3,"peak-"NR,int($score*10),".",$score,"-1",$log_FDR,int($start/2+$end/2)-$peak_start}}' \
        | bioawk -tc hdr '$5<=1000' \
        | rg -v '^#' \
        | bigtools bedtobigbed \
            -a ${template} -s start \
            - ${ref_index} fire-peaks.bb
    """
}

process PERCENT_ACCESSIBLE{
    container  'docker://tamdbho/fire:latest'
    publishDir "${params.out_dir}/${id}/trackHub/bw", mode: 'copy'
    input:
        tuple val (id), path (pileup_bed), path (pileup_tbi), val(haplotype)
        path ref_file
        path ref_index
    output:
        tuple val (id), path ("*-percent-accessible.bw")
    script:
    def suffix = params.get_hap_col_suffix(haplotype)
    def threads = 4
    """
    bgzip -cd ${pileup_bed} \
    | bioawk -tc hdr '$coverage${haplotype_suffix}>0' \
    | bioawk -tc hdr \
        'NR>1{print $1,$2,$3,100*$fire_coverage${haplotype_suffix}/$coverage${haplotype_suffix}}' \
    > tmp_${id}percent-accessible.bed

    # add fake if file is empty
    if [[ -s {output.tmp} ]]; then
        echo "File is not empty"
    else
        echo "File is empty"
        printf "${chromosome}\t0\t1\t0\\n" > {output.tmp}
    fi
    bigtools bedgraphtobigwig \
        --nzooms ${params.nzooms} -s start \
        tmp_${id}percent-accessible.bed
        ${ref_index} ${id}-${haplotype_suffix}.percent.accessible.bw
    """
}


params {
    fire = Channel.from(
    ['THP1', "${params.out_dir}/THP1/FIRE/THP1-fire-filtered.cram", "${params.out_dir}/THP1/FIRE/THP1-fire-filtered.cram.crai"],
    ['K562', "${params.out_dir}/K562/FIRE/K562-fire-filtered.cram", "${params.out_dir}/K562/FIRE/K562-fire-filtered.cram.crai"],
    ['Kasumi', "${params.out_dir}/Kasumi/FIRE/Kasumi-fire-filtered.cram", "${params.out_dir}/Kasumi/FIRE/Kasumi-fire-filtered.cram.crai"],
    ['HL60', "${params.out_dir}/HL60/FIRE/HL60-fire-filtered.cram", "${params.out_dir}/HL60/FIRE/HL60-fire-filtered.cram.crai"])



}

workflow {
// FIRE_SITES
    ref_file = file(params.reference)
    ref_index = file(params.reference_index)
    GET_CHROMS(ref_file, ref_index)
    chroms = GET_CHROMS.out
    chrom_ch = chroms.splitCsv(header: false, sep: '\n')

    // FIRE.out
    //     .combine(chrom_ch)
    //     .set { fire_sites_input_ch }

    fire = Channel.from(
    ['THP1', "${params.out_dir}/THP1/FIRE/THP1-fire-filtered.cram", "${params.out_dir}/THP1/FIRE/THP1-fire-filtered.cram.crai"],
    ['K562', "${params.out_dir}/K562/FIRE/K562-fire-filtered.cram", "${params.out_dir}/K562/FIRE/K562-fire-filtered.cram.crai"],
    ['Kasumi', "${params.out_dir}/Kasumi/FIRE/Kasumi-fire-filtered.cram", "${params.out_dir}/Kasumi/FIRE/Kasumi-fire-filtered.cram.crai"],
    ['HL60', "${params.out_dir}/HL60/FIRE/HL60-fire-filtered.cram", "${params.out_dir}/HL60/FIRE/HL60-fire-filtered.cram.crai"])

    
    fire_ch.combine(chrom_ch)
        .set { fire_sites_input_ch }

    FIRE_SITES_CHROM(fire_sites_input_ch )
    FIRE_SITES_CHROM.out
        .groupTuple()
        .set { fire_sites_ch }
    fire_sites_ch.view()
    FIRE_SITES_INDEX(FIRE_SITES.out)

// COVERAGE INFO
    BEDGRAPH(FIRE.out, ref_file, ref_index)
    COVERAGE(BEDGRAPH.out, chroms)

// FIBER LOCATIONS
    FIBER_LOCATIONS (fire_sites_input_ch)
    FIBER_LOCATIONS.out
        .groupTuple()
        .join(BEDGRAPH.out)
        .join(COVERAGE.out.minimum)
        .join(COVERAGE.out.maximum)
        .set { fiber_locations_ch }

    FIBER_LOCATIONS_INDEX(fiber_locations_ch)

// FDR
    exclude_beds = Channel.fromPath(params.bed_excludes).collect()
    shuffle_input = FIBER_LOCATIONS_INDEX.out.filtered
    SHUFFLE_EXCLUDE(shuffle_input, exclude_beds, ref_file, ref_index)
    shuffle_input
        .join(SHUFFLE_EXCLUDE.out)
        .combine(chrom_ch)
        .set { shuffle_fiber_locations_input }
    
    shuffle_fiber_locations_input.view()
    FILTERED_SHUFFLED_FIBER_LOCATIONS_CHROM(shuffle_fiber_locations_input, ref_file, ref_index)
    FILTERED_SHUFFLED_FIBER_LOCATIONS_CHROM.out.view()

    fire_sites_input_ch
        .join(FILTERED_SHUFFLED_FIBER_LOCATIONS_CHROM.out)
        .map { id, cram, crai, chrom1, chrom2, bed -> tuple(id, cram, crai, chrom1, bed) }
        .set { shuffle_pileup_chrom_input }
    shuffle_pileup_chrom_input.view()
    SHUFFLED_PILEUP_CHROM(shuffle_pileup_chrom_input)

    SHUFFLED_PILEUP_CHROM.out
        .groupTuple()
        .set { shuffle_pileups_input }
    shuffle_pileups_input.view()
    SHUFFLED_PILEUP(shuffle_pileups_input)

    SHUFFLED_PILEUP.out
        .join(COVERAGE.out.minimum)
        .join(COVERAGE.out.maximum)
        .set { fdr_table_input }   
    fdr_table_input.view()
    FDR_TABLE(fdr_table_input)

// PILEUP
    PILEUP_CHROM(fire_sites_input_ch)
    pileup_chrom_ch = PILEUP_CHROM.out
    pileup_chrom_ch.view()

    FDR_TABLE.out
        .combine(chrom_ch)
        .join(pileup_chrom_ch)
        .map { id, fdr, chrom1, chrom2, beds -> tuple(id, chrom1, beds, fdr) }
        .set { fdr_track_chrom_input }
    fdr_track_chrom_input.view()
    FDR_TRACK_CHROM(fdr_track_chrom_input)

    FDR_TRACK_CHROM.out
        .groupTuple()
        .set { pileup_input }
    pileup_input.view()
    PILEUP(pileup_input)
    pileup_beds_ch = PILEUP.out.bed
    pileup_beds_ch.view()

// FIRE PEAKS
    pileup_beds_ch
        .join(FIRE_SITES_INDEX.out)
        .combine(chrom_ch)
        .set { helper_fdr_fire_input }
    helper_fdr_fire_input.view()
    HELPER_FDR_FIRE(helper_fdr_fire_input)

// PERCENT ACCESSIBLE
    pileup_ch = PILEUP.out.bed
    pileup_ch.out
        .combine(Channel.fromList(params.not_unk))
        .set { percent_accessible_input }
    PERCENT_ACCESSIBLE(percent_accessible_input, ref_file, ref_index)
}