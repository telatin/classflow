process PIGZ_COMPRESS {
    tag "${sample_id}@${task.cpus}"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::seqfu=1.17.0' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.17.0--hbd632db_0' :
        'quay.io/biocontainers/seqfu:1.17.0--hbd632db_0' }"
    
    input:
    tuple val(sample_id), path(file) 

    
    output:
    tuple val(sample_id), path("*gz")  

    script:
    """
    cat ${file} | pigz -p ${task.cpus} > ${sample_id}.gz
    """     
}

process INTERLEAVE {
    tag "${sample_id}@${task.cpus}"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::seqfu=1.17.0' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.17.0--hbd632db_0' :
        'quay.io/biocontainers/seqfu:1.17.0--hbd632db_0' }"
    
    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("${sample_id}.fq.gz")  

    script:
    """
    seqfu interleave -1 ${reads[0]} -2 ${reads[1]} | gzip -c > ${sample_id}.fq.gz
    """ 
}
process DEINTERLEAVE {
    tag "${sample_id}@${task.cpus}"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::seqfu=1.17.0' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.17.0--hbd632db_0' :
        'quay.io/biocontainers/seqfu:1.17.0--hbd632db_0' }"


    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.fq.gz")  

    script:
    """
    seqfu deinterleave -o ${sample_id} ${reads} 
    gzip ${sample_id}_R{1,2}.fq
    """ 
}
