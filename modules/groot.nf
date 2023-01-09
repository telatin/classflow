
process GROOT {
    //groot:1.1.2--hef68116_1  
    conda (params.enable_conda ? 'bioconda::groot=1.1.2' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/groot:1.1.2--hef68116_1' :
        'quay.io/biocontainers/groot:1.1.2--hef68116_1' }"
    
    label 'process_medium'
    tag "${sample_id}@${task.cpus}"


    publishDir "$params.outdir/amr/", 
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    path(grootdb)
   
    
    output:
    tuple val(sample_id), path("*.txt")

 
    script:
    """  
    groot version > groot_version.log
    gzip -dc "${reads[0]}" | groot align -i "$grootdb" -p  24 | groot report -c  0.97 > ${sample_id}.txt
    """
}