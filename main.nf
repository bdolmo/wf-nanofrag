#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameters for input/output directories and paths
params.tumor = null                                        // Path to tumor BAM list (txt file)
params.normal = null                                       // Path to normal BAM list (txt file)
params.reference = null                                         // Path to the reference genome (FASTA file)
params.out_dir = "./nanofrag_results"                          // Directory to save output metrics
params.threads = 12                                             // Number of threads
params.skip_small_variants = true                              // Whether to skip small variants (flag)

if (params.help) {
    log.info """
    Nanofrag Analysis
    =============================================
    Usage:
        nextflow run main.nf --tumor <file> --normal <file> --reference <file> --out_dir <dir> --threads <number>
    Input:
        * --tumor: directory with BAM files
        * --normal: directory with normal BAM files
        * --reference: path to the reference genome (FASTA file)
        * --out_dir: output directory. Default [${params.out_dir}]
        * --threads: number of threads. Default [${params.threads}]
        * --skip_small_variants: flag to skip small variants. Default [${params.skip_small_variants}]
    """
    exit 0
}

// Ensure mandatory parameters are provided
if (!params.tumor || !params.normal || !params.reference) {
    throw new IllegalArgumentException("Missing required parameters: --tumor, --normal, or --reference")
}

// Ensure the output directory exists
if (!new File(params.out_dir).exists()) {
    new File(params.out_dir).mkdirs()
}

// Create channels for inputs
tumor_list = Channel.fromPath(params.tumor)
normal_list = Channel.fromPath(params.normal)
reference = Channel.fromPath(params.reference)



def tumor_input = file(params.tumor, stageAs:"tumor_path")
def normal_input = file(params.normal, stageAs:"normals_path")
def ref_input = file(params.reference, stageAs:"reference_path")


def tumor_basename = file(params.tumor).getBaseName()
def normal_basename = file(params.normal)getBaseName()
def ref_basename = file(params.reference).getBaseName()

def tumor_dirname = file(params.tumor).parent 
def normal_dirname = file(params.normal).parent 
def ref_dirname = file(params.reference).parent 


def nanofrag_dir = file(params.nanofrag).parent


process runNanofrag {
    tag "nanofrag"
    memory '12 GB'
    cpus params.threads

    input:
    // path tumor_list
    // path normal_list
    // path reference

    output:

    script:
    """
    echo "Running Nanofrag using docker run command..."
    docker run  -i --rm \\
        -v ${nanofrag_dir}/:/nanofrag_script/ \\
        -v ${normal_input}/:/tumors/ \\
        -v ${params.normal}/:/normals/ \\
        -v ${params.reference}/:/ref_dir/\\
        -v ${params.out_dir}:/out_dir/ \\
        bdolmo/python_env_nanofrag:latest /nanofrag_script/nanofrag.py \\
        --docker_output ${params.out_dir} \\
        --tumor_list /tumors/ \\
        --normal_list /normals/ \\
        --reference /ref_dir/ \\
        --output_dir /out_dir/ \\
        --threads ${task.cpus} \\
        ${params.skip_small_variants ? '--skip_small_variants' : ''}
    """
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
workflow {
    runNanofrag()
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
