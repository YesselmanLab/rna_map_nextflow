#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// TODO make sure num of chunks is read into a program somewhere

process split_fastqs {
    input:
        tuple val(name), path(r1_path), path(r2_path)

    output:
        path("*.fastq.gz")

    script:
    """
    rna-map-tools split-fastqs $r1_path $r2_path $params.split_fastq_chunks
    """
}


process trim_galore {
    input:
        tuple file(r1_chunk), file(r2_chunk)

    output:
        tuple file("*R1*.fq.gz"), file("*R2*.fq.gz") 

    script:
    """
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
    rna-map-tools demultiplex ${params.data_csv} $fastq_chunk_1 $fastq_chunk_2
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
    // check to make sure all required parameters are specified
    if (params.fastq_dirs.length() == 0) {
        println("No fastq directories specified")
        exit 1
    }
    if (params.data_csv.length() == 0) {
        println("No data csv specified")
        exit 1
    }
    l = params.fastq_dirs.tokenize(',')
    l = l.collect{ it + "/*{1,2}.fastq.gz" }
    fastq_dir_ch = Channel.fromFilePairs(l, flat: true)
    //fastq_dir_ch = fastq_dir_ch.map { it -> it + "/*.fastq.gz" }
    //fastq_dir_ch = Channel.fromFilePairs(fastq_dir_ch, flat: true)
    //fastq_dir_ch.view()
    // load fastq files from directories 
    // should we split? 
    if (params.split_fastqs) {
        // spit fastqs into chunks
        fastq_dir_ch = fastq_dir_ch | split_fastqs
        fastq_dir_ch = fastq_dir_ch.flatten()
        fastq_dir_ch = fastq_dir_ch.map{ it-> tuple(it.toString().tokenize("_")[3], it)}
        fastq_dir_ch = fastq_dir_ch.groupTuple(by: 0)
        fastq_dir_ch.view()

    }
    /*else {
        println("made it")
        fastq_dir_ch = Channel.fromFilePairs('fastqs/*.fastq.gz', flat: true)
        fastq_dir_ch.view()
    }*/   
    /*
    // run trim galore on each chunk
    fastq_dir_ch.view()
    trim_fastqs_ch = fastq_dir_ch | trim_galore
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
    */
}  
