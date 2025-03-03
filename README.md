# nf-fiberseq

Nextflow workflow running ft predict m6a, pbsv, deepvariant, sniffles, hiphase, fire, pb-CpG-calling on hifi bam files. 

Adapted from Stergachis lab's snakemake workflow

(https://github.com/mrvollger/k-mer-variant-phasing/blob/master/workflow/rules/kmer.smk

https://github.com/furlan-lab/FIRE/blob/main/workflow/rules/fire-peaks.smk)

## Workflow

1. Run ft predict-m6a
2. Align using pbmm2 after filtering using zmwfilter (PacBio). 
3. Extract _m6a, _msp, _nuc, and _cpg bed files
4. Sort and index BAM files. Optional: merge BAM files with same sample ID
5. Variant calling using: deepvariant (small variants), pbsv and sniffles (structural variants)
6. Phase variant calls using HiPhase
7. Run ft fire using haplotagged BAM 
8. Run align_bam_to_cpg_scores using haplotagged BAM

# Components

## Input data

CSV file listing sample ID and corresponding paths to all input BAM files

Example:
```
id,reads_bam
sampleA,/path/to/hifi/sampleA.bam
sampleB,/path/to/hifi/sampleB.bam
```
In this directory, there is a make_samplesheet.py that might be helpful to generate samplesheet.csv if sample id is not included in file name

## Nextflow.config

Find nextflow.congig file and specify parameters that will be used in nextflow pipeline:
- ```sample_sheet``` = '/path/to/samplesheet.csv'
- ```out_dir``` = '/path/to/output/directory' (important output will be copied here from workDir)
- ```reference``` = '/path/to/reference.fasta'
- ```reference_index``` = 'path/to/reference.fasta.fai'
- ```temp_dir``` = '/path/to/temp/directory' (this directory will not be removed upon workflow completion , but exists as symlink from workDir)

- ```bam_merge``` = true/false (set true if there are multiple BAM files per sample and need merging using samtools merge)
- ```workDir``` = Nextflow built in cache system, specify directory where all the outputs from various runs can be stored in

Other parameters are most likely fixed - defaults from Stergachis lab workflow. 

## Containers

Apptainer is enabled for this nextflow project. Software used in each process is either outsourced from public container or built from attached Dockerfiles (found in docker folder)

# Running the workflow

Attached is NF-main.sh script which includes:

```
nextflow \
    run main.nf \
    -c ${NXF_CONFIG} -profile slurm \
    -with-report nextflow.report.html \
    -resume
```
This script automatically runs main.nf file which contains the main workflow and send each process to HPC SLURM. 
```-profile``` argument can be changed to 'local' for local execution. Otherwise, sbatch this script as is to keep default configurations. 

# Output files
In the specified output directory expect 1 folder for each sample ID which should include:

```
m6a/            # m6a BAM files (unaligned)
BAM/            # sorted/merged m6A BAM (unaligned) and haplotagged BAM files
deepvariant/    # gvcf.gz and vcf.gz output files (indexed) from running deepvariant
FIRE/           # cram and crai files from running ft fire
ft_extract/     # bed files extracted from m6a BAM
hiphase/        # deepvariant and pbsv phased vcf files 
pbmm2/          # aligned BAM files
pbsv/           # vcf.gz output files (indexed) from running pbsv
sniffles/       # vcf file from running sniffles
```

Work in progress: visual.nf to visualize FIRE outputs

# Authors

Analysis code was written by Stergachis lab. Workflow code was written by Tam Ho (tam.db.ho@gmail.com)
