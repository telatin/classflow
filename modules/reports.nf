process MULTIQC {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::multiqc=1.12" : null)   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    publishDir "$params.outdir/", 
        mode: 'copy'
        
    input:
    path("*")

    output:
    path("report.html")

    script:
    """
    multiqc --force --filename report.html --no-data-dir --title "Classflow" --comment "Summary of findings by ClassFlow" .
    """

}