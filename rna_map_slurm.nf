#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastq_dirs = []
params.barcode_csv = "$baseDir/inputs/data.csv"
params.split_chunks = 200

process split_fastqs {
    input:
        path(fastq_dir)

    output:
        file("*.fastq.gz")

    script:
    """
    source $baseDir/setup_env
    python $baseDir/scripts/split_fastqs.py $fastq_dir $params.split_chunks
    """
}

process trim_galore {
    input:
        tuple file(r1_chunk), file(r2_chunk)

    output:
        tuple file("*R1*.fq.gz"), file("*R2*.fq.gz") 

    script:
    """
    source $baseDir/setup_env
    trim_galore --fastqc --quality 0 --paired $r1_chunk $r2_chunk
    """
}

process demultiplex_fastq {
    input:
        tuple path(fastq_chunk_1), path(fastq_chunk_2)

    output:
        path "[ACGT]*"

    script:
    """
    source $baseDir/setup_env
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
    source $baseDir/setup_env
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
    source $baseDir/setup_env
    python ${baseDir}/scripts/join_zip_files_multi.py ${barcode} ${barcode_path}
    """
}

process run_rna_map {
    input:
    path(barcode_dir)
    
    output:
    path("output-*")

    script:
    """
    source $baseDir/setup_env
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
    source $baseDir/setup_env
    python ${baseDir}/scripts/combine_output_final.py ${params.barcode_csv} ${barcode} ${rna_map_dirs}
    """

}

workflow {
    // load fastq files from directories 
    fastq_dir_ch = Channel.fromPath(params.fastq_dirs.tokenize(','))
    // split fastq files into chunks
    split_fastqs_ch = fastq_dir_ch | split_fastqs
    split_fastqs_ch = split_fastqs_ch.flatten()
    r1_ch = split_fastqs_ch.filter { it.toString().contains('R1') }
    r2_ch = split_fastqs_ch.filter { it.toString().contains('R2') }
    split_fastqs_ch = r1_ch.merge(r2_ch)
    // run trim galore on each chunk
    trim_fastqs_ch = split_fastqs_ch | trim_galore
    // demultiplex each chunk with sabre 
    barcode_file_ch = trim_fastqs_ch | demultiplex_fastq
    barcode_file_ch = barcode_file_ch.flatten()
    // demultiplex each chunk with internal script
    int_demult_ch = barcode_file_ch | internal_demultiplex
    grouped_demultiplexed_ch = int_demult_ch.
        map({file -> 
            def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
        }).
        groupTuple(by: 0)
    // join zip files 
    joined_dirs = grouped_demultiplexed_ch | join_zip_files
    joined_dirs = joined_dirs.flatten()
    // run rna map on each construct individually
    rna_map_outputs_ch = joined_dirs | run_rna_map
    // combine output files and generate final output
    grouped_rna_outputs_ch = rna_map_outputs_ch.
        map({file -> 
            def key = file.name.toString().tokenize('-').get(1)
            return tuple(key, file)
        }).
        groupTuple(by: 0)
    grouped_rna_outputs_ch | combine_output_final
}  
