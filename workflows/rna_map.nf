#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process split_fastqs {
    input:
        tuple val(name), path(r1_path), path(r2_path)

    output:
        path("*.fastq.gz")

    script:
    """
    rna-map-nextflow split-fastqs $r1_path $r2_path $params.split_fastq_chunks
    """
}


// runs trim galore in a temporary directory and moves the outputed fastq files to the current directory
// this is to reduce the number of files created in $WORK on the cluster
process trim_galore {
    input:
        tuple val(name), path(r1_chunk), path(r2_chunk)

    output:
        tuple file("*R1*.fq.gz"), file("*R2*.fq.gz") 

    script:
    """
    random_string=\$(cat /dev/urandom | env LC_ALL=C tr -dc 'A-Za-z0-9' | head -c 10)
    outdir="$params.tmp_dir/\$random_string"
    mkdir -p \$outdir
    trim_galore --fastqc --quality 24 --length 60 --paired -o \$outdir  $r1_chunk $r2_chunk 
    mv \${outdir}/*R1*.fq.gz .
    mv \${outdir}/*R2*.fq.gz .
    rm -rf \$outdir
    """
}

process demultiplex_fastq {
    input:
        tuple path(fastq_chunk_1), path(fastq_chunk_2)

    output:
        path "[ACGT]*"

    script:
    """
    rna-map-nextflow demultiplex $params.data_csv $fastq_chunk_1 $fastq_chunk_2
    """
}

process internal_demultiplex {
    input:
        path(barcode_dir)
    
    output:
        path("*.demultiplexed.zip")
    
    script:
    """
    rna-map-nextflow int-demultiplex $params.input_dir $barcode_dir -o $params.tmp_dir
    """
}

process join_zip_files {
    input:
    tuple val(barcode), path(barcode_path)
    
    output:
    path("output/[ACGT]*")

    script:
    """
    rna-map-nextflow join-int-demultiplex-zips --threads $task.cpus $barcode $barcode_path -o $params.tmp_dir
    """
}

process run_rna_map {
    input:
    path(barcode_dirs)
    
    output:
    path("output-*")

    script:
    """
    for barcode_dir in $barcode_dirs; do
        rna-map-nextflow run-rna-map $params.input_dir \$barcode_dir -o $params.tmp_dir
    done 
    """
}

process combine_output_final {
    publishDir path: "${launchDir}/processed", mode:'copy', overwrite: true

    input: 
    tuple val(barcode), path(rna_map_dirs)

    output:
    path("*")

    script:
    """
    rna-map-nextflow combine-rna-map-outputs $barcode $params.input_dir $rna_map_dirs
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
    // TODO would be nice to check to make sure there are actually fastqs? 
    l = params.fastq_dirs.tokenize(',')
    l = l.collect{ it + "/*{1,2}.fastq.gz" }
    fastq_dir_ch = Channel.fromFilePairs(l, flat: true)
    // load fastq files from directories 
    // should we split? 
    if (params.split_fastqs) {
        // spit fastqs into chunks
        fastq_dir_ch = fastq_dir_ch | split_fastqs
        fastq_dir_ch = fastq_dir_ch.flatten()
        fastq_dir_ch = fastq_dir_ch.map{ it-> tuple(it.toString().tokenize("_")[-2], it)}
        fastq_dir_ch = fastq_dir_ch.groupTuple(by: 0)
        fastq_dir_ch = fastq_dir_ch.map{ it-> tuple(it[0], it[1][0], it[1][1])}
    }   
    // run trim galore on each chunk
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
    joined_dirs = joined_dirs.collate(4, false)
    // run rna map on each construct individually
    rna_map_outputs_ch = joined_dirs | run_rna_map
    rna_map_outputs_ch = rna_map_outputs_ch.flatten()
    // combine output files and generate final output
    grouped_rna_outputs_ch = rna_map_outputs_ch.
        map({file -> 
            def key = file.name.toString().tokenize('-').get(1)
            return tuple(key, file)
        }).
        groupTuple(by: 0)
    grouped_rna_outputs_ch | combine_output_final
}  
