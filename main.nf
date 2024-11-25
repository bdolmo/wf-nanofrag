#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameters for input/output directories and paths
params.tumor_list = null                                        // Path to tumor BAM list (txt file)
params.normal_list = null                                       // Path to normal BAM list (txt file)
params.reference = null                                         // Path to the reference genome (FASTA file)
params.out_dir = "./nanofrag_results"                          // Directory to save output metrics
params.threads = 12                                             // Number of threads
params.skip_small_variants = true                              // Whether to skip small variants (flag)

if (params.help) {
    log.info """
    Nanofrag Analysis
    =============================================
    Usage:
        nextflow run main.nf --tumor_list <file> --normal_list <file> --reference <file> --out_dir <dir> --threads <number>
    Input:
        * --tumor_list: path to the text file listing tumor BAM files
        * --normal_list: path to the text file listing normal BAM files
        * --reference: path to the reference genome (FASTA file)
        * --out_dir: output directory. Default [${params.out_dir}]
        * --threads: number of threads. Default [${params.threads}]
        * --skip_small_variants: flag to skip small variants. Default [${params.skip_small_variants}]
    """
    exit 0
}

// Ensure mandatory parameters are provided
if (!params.tumor_list || !params.normal_list || !params.reference) {
    throw new IllegalArgumentException("Missing required parameters: --tumor_list, --normal_list, or --reference")
}

// Ensure the output directory exists
if (!new File(params.out_dir).exists()) {
    new File(params.out_dir).mkdirs()
}

// Create channels for inputs
tumor_list = Channel.fromPath(params.tumor_list)
normal_list = Channel.fromPath(params.normal_list)
reference = Channel.fromPath(params.reference)

process runNanofrag {
    tag "nanofrag"
    memory '12 GB'
    cpus params.threads
    conda params.conda_yaml

    input:
    path tumor_list
    path normal_list
    path reference

    output:

    script:
    """
    echo "Running Nanofrag with tumor list ${tumor_list}, normal list ${normal_list}, and reference genome ${reference}..."
    python3 ${projectDir}/bin/nanofrag/nanofrag.py \\
        --tumor_list ${tumor_list} \\
        --normal_list ${normal_list} \\
        --reference ${reference} \\
        --output_dir . \\
        --threads ${params.threads} \\
        ${params.skip_small_variants ? '--skip_small_variants' : ''}
    """
}

// Workflow definition
workflow {
    runNanofrag(tumor_list, normal_list, reference)
}

workflow.onComplete {
    println "Nanofrag analysis completed successfully. Results are saved in '${params.out_dir}'."
}




// process fragmentExtractor {
//     tag "${bam_file.simpleName}"                                // Tag the process with the BAM filename
//     memory '2 GB'

//     input:
//     path bed_file
//     path bam_file
//     path bai_file

//     output:
//     path "${bam_file.simpleName}_fragments.bed" // Output BED file
//     path "${bam_file.simpleName}_fragments.wig" // Output WIG file

//     script:
//     """
//     echo "Processing ${bam_file} with ${bed_file}..."
//     fragment_extractor ${bam_file} ${bed_file} \
//         ${bam_file.simpleName}_fragments.bed \
//         ${bam_file.simpleName}_fragments.wig
//     """
// }

// // Workflow definition
// workflow {
//     results = fragmentExtractor(windows_bed, bam_files, bai_files)

// }

// workflow.onComplete {
//     println "Workflow execution completed successfully. Fragmentation metrics are saved in the '${params.out_dir}' directory."
// }
