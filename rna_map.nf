#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_dir = "./fastq_input"
params.output_dir = "./demultiplexed_output"
params.barcodes = "./barcodes.txt"
params.barcode_csv = "$baseDir/data.csv"

// need to make split fastq a job
// need to run trim galore first
// need to combine runs at the end

process demultiplex_fastq {
    input:
        tuple val(filename), path(fastq_chunk_1), path(fastq_chunk_2)

    output:
        path "[ACGT]*"

    script:
    """
    python ${baseDir}/scripts/sabre_demultiplex.py ${params.barcode_csv} ${fastq_chunk_1} ${fastq_chunk_2}
    """
}

process internal_demultiplex {
    input:
        path(barcode_dir)
    
    output:
        path("*.demultiplexed.zip")
    
    script:
    """
    python ${baseDir}/scripts/internal_demultiplex.py ${params.barcode_csv} ${barcode_dir}
    """
}

process join_zip_files {
    input:
    tuple val(barcode), path(barcode_path)
    
    output:
    path("output/[ACGT]*")

    script:
    """
    python ${baseDir}/scripts/join_zip_files.py ${barcode} ${barcode_path}
    """
}

process run_rna_map {
    input:
    path(barcode_dir)
    
    output:
    path("output-*")

    script:
    """
    python ${baseDir}/scripts/run_rna_map.py ${baseDir} ${barcode_dir}
    """
}

process combine_output_final {
    publishDir path: "${baseDir}/processed", mode:'copy', overwrite: true

    input: 
    tuple val(barcode), path(rna_map_dirs)

    output:
    path("*")

    script:
    """
    python ${baseDir}/scripts/combine_output_final.py ${params.barcode_csv} ${barcode} ${rna_map_dirs}
    """

}

workflow {
    split_fasta_ch = \
        channel.fromFilePairs("$baseDir/test_R{1,2}.fastq.gz", flat: true) | \
        splitFastq(by: 10_000, pe: true, file: true, compress: true)
    barcode_file_ch = split_fasta_ch | demultiplex_fastq
    barcode_file_ch = barcode_file_ch.flatten()
    int_demult_ch = barcode_file_ch | internal_demultiplex
    grouped_demultiplexed_ch = int_demult_ch.
        map({file -> 
            def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
        }).
        groupTuple(by: 0)
    joined_dirs = grouped_demultiplexed_ch | join_zip_files
    joined_dirs = joined_dirs.flatten()
    rna_map_outputs_ch = joined_dirs | run_rna_map
    grouped_rna_outputs_ch = rna_map_outputs_ch.
        map({file -> 
            def key = file.name.toString().tokenize('-').get(1)
            return tuple(key, file)
        }).
        groupTuple(by: 0)
    grouped_rna_outputs_ch | combine_output_final
}  
