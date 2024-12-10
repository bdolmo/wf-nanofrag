#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameters for input/output directories and paths
params.tumor = null                                        // Path to tumor BAM directory
params.normal = null                                       // Path to normal BAM directory
params.reference = null                                    // Path to the reference genome (FASTA file)
params.out_dir = "./nanofrag_results"                     // Directory to save output metrics
params.threads = 4                                        // Number of threads
params.skip_small_variants = true                         // Whether to skip small variants (flag)

if (params.help) {
    log.info """
    Nanofrag Analysis
    =============================================
    Usage:
        nextflow run main.nf --tumor <directory> --normal <directory> --reference <file> --out_dir <dir> --threads <number>
    Input:
        * --tumor: directory with tumor BAM files
        * --normal: directory with normal BAM files
        * --reference: path to the reference genome (FASTA file)
        * --out_dir: output directory. Default [${params.out_dir}]
        * --threads: number of threads. Default [${params.threads}]
        * --skip_small_variants: flag to skip small variants. Default [${params.skip_small_variants}]
    """
    exit 0
}

// Ensure mandatory parameters are provided
if (!params.tumor || !params.normal || !params.reference || !params.out_dir) {
    throw new IllegalArgumentException("Missing required parameters: --tumor, --normal, --reference, or --out_dir ")
}

// Ensure the output directory exists
if (!new File(params.out_dir).exists()) {
    new File(params.out_dir).mkdirs()
}

// Create input channels
tumor_dir = params.tumor
normal_dir = params.normal
reference_file = params.reference
reference_fai = params.reference + ".fai"
output_dir = params.out_dir

println(reference_file)

process runNanofrag {
    tag "nanofrag"
    memory '15 GB'
    cpus params.threads
    container 'bdolmo/nanofrag:latest'  // Use Nextflow Docker integration

    input:
    path tumor_dir
    path normal_dir
    path reference_file
    path reference_fai
    path output_dir


    output:

    script:
    """
    python3.12 /opt/nanofrag/nanofrag.py \\
        --docker_output ${output_dir} \\
        --tumor_list ${tumor_dir} \\
        --normal_list ${normal_dir} \\
        --reference ${reference_file} \\
        --output_dir ${output_dir} \\
        --threads ${task.cpus} \\
        ${params.skip_small_variants ? '--skip_small_variants' : ''}
    """
}

// Workflow definition
workflow {
    runNanofrag(tumor_dir, normal_dir, reference_file, reference_fai, output_dir)
}

workflow.onComplete {
    println "Nanofrag analysis completed successfully. Results are saved in '${params.out_dir}'."
}



// process runNanofrag {
//     tag "nanofrag"
//     memory '12 GB'
//     cpus params.threads
//     docker bdolmo/python_env_nanofrag

//     input:
//     path tumor_list
//     path normal_list
//     path reference

//     output:

//     script:
//     """
//     # conda env create -f ${params.conda_yaml}
//     # source activate base
//     # conda init
//     # conda activate nanofrag_env
//     echo "Running Nanofrag with tumor list ${tumor_list}, normal list ${normal_list}, and reference genome ${reference}..."
//     nanofrag.py \\
//         --tumor_list ${tumor_list} \\
//         --normal_list ${normal_list} \\
//         --reference ${reference} \\
//         --output_dir . \\
//         --threads ${params.threads} \\
//         ${params.skip_small_variants ? '--skip_small_variants' : ''}
//     """
// }

// Workflow definition


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

