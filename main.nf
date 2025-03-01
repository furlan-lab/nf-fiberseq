#!/usr/bin/env nextflow

log.info """\
FIBERSEQ
==========================================
samplesheet      : $params.sample_sheet
output_directory : $params.out_dir
temp_directory   : $params.temp_dir
ref_name         : $params.ref_name

FIRE
=========================================
filter_flags        : $params.filter_flag
min_msp             : $params.min_msp
min_ave_msp_size    : $params.min_ave_msp_size
min_contig_length   : $params.min_contig_length
coverage_within_n_sd: $params.coverage_within_n_sd
min_coverage        : $params.min_coverage
min_per_acc_peak    : $params.min_per_acc_peak
min_frac_accessible : $params.min_frac_accessible
max_frac_overlap    : $params.max_frac_overlap
max_peak_fdr        : $params.max_peak_fdr
min_fire_fdr        : $params.min_fire_fdr
"""

// Import processes
include { 
    FT_M6A;
    PBI_INDEX;
    ZMW_FILTER;
    ZMW_INCLUDE;
    PBMM2_ALIGN;
    FT_EXTRACT;
    BED_SORT;
    BAM_SORT;
    BAM_MERGE;
    BAM_INDEX;
    DEEPVARIANT;
    DEEPVARIANT_INDEX;
    PBSV_DISCOVER;
    PBSV;
    PBSV_INDEX;
    SNIFFLES;
    HIPHASE;
    FIRE;
    pb_CpG_calling
} from './modules.nf'

include { BAM_INDEX as BAM_INDEX2 } from './modules.nf' 

workflow {
    Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.id, file(row.reads_bam)) } 
    .set { bam_files }

    ref_file  = file(params.reference)
    ref_index = file(params.reference_index)
// ft predict-m6a
    FT_M6A(bam_files)
    PBI_INDEX(FT_M6A.out)
// align
    ZMW_FILTER(PBI_INDEX.out)
    ZMW_INCLUDE(
        PBI_INDEX.out, 
        ZMW_FILTER.out)
    PBMM2_ALIGN(
        ZMW_INCLUDE.out,
        ref_file,
	    ref_index)
// bed files
    FT_EXTRACT(PBMM2_ALIGN.out)

    bed_ch = FT_EXTRACT.out
                        .flatMap { id, paths -> 
                        paths.collect { path -> tuple(id, path) }}
    bed_ch.view()
    BED_SORT(bed_ch)
    // BED_INDEX(BED_SORT.out)

    BAM_SORT(PBMM2_ALIGN.out)
    sortedBAM_ch = BAM_SORT.out.groupTuple()
    sortedBAM_ch.view()
    
    // If BAM files need to be merged (one sample = multiple BAM files) edit params.bam_merge to true (nextflow.config)
    if (params.bam_merge) {
        BAM_MERGE(sortedBAM_ch)
        BAM_INDEX(BAM_MERGE.out)
    } else {
        BAM_INDEX(BAM_SORT.out)
    }
    
// DEEPVARIANT
    DEEPVARIANT(BAM_INDEX.out, 
		ref_file, 
		ref_index)
    DEEPVARIANT_INDEX(DEEPVARIANT.out)

// PBSV
    // Identify signatures of structural variation. 
    // Reduces all aligned reads to those that are relevant to calling structural variants. 
    PBSV_DISCOVER(BAM_INDEX.out, ref_file,ref_index)
        .set { pbsv_ch }

    BAM_INDEX.out
        .join(pbsv_ch)
        .set { pbsv_input_ch }
    pbsv_input_ch.view()

    PBSV(pbsv_input_ch, 
        ref_file, 
        ref_index)
    PBSV_INDEX(PBSV.out)
    PBSV_INDEX.out.view()

// Sniffles
    SNIFFLES(BAM_INDEX.out, 
            ref_file, 
            ref_index)

// HiPhase
    deepvariant_vcf = DEEPVARIANT_INDEX.out.vcf
    pbsv_output = PBSV_INDEX.out
    BAM_INDEX.out
        .join(deepvariant_vcf)
        .join(PBSV_INDEX.out)
        .set { hiphase_input_ch }
    hiphase_input_ch.view()
   
    HIPHASE(hiphase_input_ch,
            ref_file,
            ref_index)
// FIRE
    haplotagged_bam = HIPHASE.out.bam
    haplotagged_bam.view()
    FIRE(haplotagged_bam, ref_file)
    FIRE.out.view()

// pb_CpG_calling
    BAM_INDEX2(haplotagged_bam)
    BAM_INDEX2.out
        .set{pb_CpG_calling_input}
    pb_CpG_calling(pb_CpG_calling_input, ref_file, ref_index)
}

