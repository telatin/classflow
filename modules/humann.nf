process CHECK_MPA {

  conda (params.enable_conda ? 'bioconda::humann=3.6' : null)    
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann%3A3.6--pyh7cba7a3_1' :
        'quay.io/biocontainers/humann%3A3.6--pyh7cba7a3_1' }"
  
  publishDir "$params.outdir/humann/", 
        mode: 'copy'

  input:
  path(metaphlandb)
  val(mpa_ver)
  
  output:
  path("mpa_version.txt")

  script:
  """
  MPA=\$(cat ${metaphlandb}/mpa_latest)
  if [ "\${MPA}" != "${mpa_ver}" ]; then
    echo "MPA version mismatch: ${mpa_ver} != \${MPA}" > mpa_version.txt
    exit 1
  else
    echo "\${MPA}" > mpa_version.txt
  fi
  """
}
process HUMANN {
    conda (params.enable_conda ? 'bioconda::humann=3.6' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann%3A3.6--pyh7cba7a3_1' :
        'quay.io/biocontainers/humann:3.6--pyh7cba7a3_1' }"
    
    tag "${sample_id}@${task.cpus}"
    label 'humann'

    publishDir "$params.outdir/humann/", 
        mode: 'copy',
        saveAs: { file -> file.take(file.lastIndexOf('.')) }

    input:
    tuple val(sample_id), path(reads) 
    path(chocophlan)
    path(uniref)
    path(metaphlandb)
    val(mpa_ver)
    path("checked.txt")
    
    output:
    tuple val(sample_id), path("*.tsv")

 
    script:
    /*
    usage: metaphlan --input_type {fastq,fasta,bowtie2out,sam} [--force]
                   [--bowtie2db METAPHLAN_BOWTIE2_DB] [-x INDEX]
                   [--bt2_ps BowTie2 presets] [--bowtie2_exe BOWTIE2_EXE]
                   [--bowtie2_build BOWTIE2_BUILD] [--bowtie2out FILE_NAME]
                   [--min_mapq_val MIN_MAPQ_VAL] [--no_map] [--tmp_dir]
                   [--tax_lev TAXONOMIC_LEVEL] [--min_cu_len]
                   [--min_alignment_len] [--add_viruses] [--ignore_eukaryotes]
                   [--ignore_bacteria] [--ignore_archaea] [--ignore_ksgbs]
                   [--ignore_usgbs] [--stat_q] [--perc_nonzero]
                   [--ignore_markers IGNORE_MARKERS] [--avoid_disqm] [--stat]
                   [-t ANALYSIS TYPE] [--nreads NUMBER_OF_READS]
                   [--pres_th PRESENCE_THRESHOLD] [--clade] [--min_ab]
                   [-o output file] [--sample_id_key name]
                   [--use_group_representative] [--sample_id value]
                   [-s sam_output_file] [--legacy-output] [--CAMI_format_output]
                   [--unclassified_estimation] [--mpa3] [--biom biom_output]
                   [--mdelim mdelim] [--nproc N] [--subsampling SUBSAMPLING]
                   [--subsampling_seed SUBSAMPLING_SEED] [--install] [--offline]
                   [--force_download] [--read_min_len READ_MIN_LEN] [-v] [-h]
                   [INPUT_FILE] [OUTPUT_FILE]
    */
    """  
    humann --version 2>&1 > humann.version
    humann \\
      -i "${reads}" -o "out" \\
      --nucleotide-database ${chocophlan} \\
      --protein-database ${uniref} \\
      --metaphlan-options=" --index ${mpa_ver} --bowtie2db ${metaphlandb} " \\
      --threads ${task.cpus}
    mv out/*_humann_temp/*_metaphlan_bugs_list.tsv out/
    rm -rf out/*_humann_temp
    mv out/* .
    """
}
process HUMANN_MERGE {
    conda (params.enable_conda ? 'bioconda::humann=3.6' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann%3A3.6--pyh7cba7a3_1' :
        'quay.io/biocontainers/humann:3.6--pyh7cba7a3_1' }"
    
    tag "${sample_id}@${task.cpus}"
    label 'humann'

    publishDir "$params.outdir/tables/", 
        mode: 'copy'
        

    input:
    path("*")

    
    output:
    path("metaphlan.tsv")
    path("gene_families.tsv")
    path("pathway_abundance.tsv")
 
    script:
    """  
    merge_metaphlan_tables.py *bugs_list.tsv > metaphlan.tsv
    humann_join_tables -i . -o gene_families.tsv     --file_name  _genefamilies
    humann_join_tables -i . -o pathway_abundance.tsv --file_name  _pathabundance
    """
}

process HUMANN_LEGACY {
    conda (params.enable_conda ? 'bioconda::humann=3.6' : null)    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.6--pyh7cba7a3_1' :
        'quay.io/biocontainers/humann:3.6--pyh7cba7a3_1' }"
    
    tag "${sample_id}@${task.cpus}"
    label 'humann'

    publishDir "$params.outdir/humann/", 
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    path(chocophlan)
    path(uniref)
    path(metaphlandb)
    val(mpa_ver)
    
    output:
    path("out")

 
    script:
    """
    humann --version 2>&1 > humann.version
    humann \\
      -i "${reads}" -o "out" \\
      --nucleotide-database ${chocophlan} \\
      --protein-database ${uniref} \\
      --metaphlan-options=" --index ${mpa_ver} --bowtie2db ${metaphlandb} --unknown_estimation " \\
      --threads ${task.cpus}
    """  
} 