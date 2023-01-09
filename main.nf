/* 
 *   Input parameters 
 */
nextflow.enable.dsl = 2
params.dir = "$baseDir/nano/"
params.pattern = "*_R{1,2}.fastq.gz"
def reads = params.dir + params.pattern
 
params.outdir = "$baseDir/denovo"

params.enable_conda = false
params.krakendb    = false
params.chocophlan  = false
params.uniref      = false
params.metaphlandb = false
params.mpa_ver     = false      
params.grootdb     = false
// prints to the screen and to the log
log.info """
         GMH Classification (version 1)
         ===================================
         dir,pattern  : ${params.dir}${params.pattern}
         outdir       : ${params.outdir}
         
         krakendb     : ${params.krakendb}
         metaphlandb  : ${params.metaphlandb}
         uniref       : ${params.uniref}
         chocophlan   : ${params.chocophlan}
         groot        : ${params.grootdb}
         """
         .stripIndent()

/* 
   check reference path exists 
*/
 
def krakenDbPath = file(params.krakendb, checkIfExists: true)
def unirefDB     = file(params.uniref, checkIfExists: true)
def chocophlanDB = file(params.chocophlan, checkIfExists: true)
def metaphlanDB  = file(params.metaphlandb, checkIfExists: true)
file("${params.krakendb}/hash.k2d", checkIfExists: true)

/*
  Modules
*/

include { KRAKEN2_HOST; KRAKEN2_REPORT }                from './modules/kraken'
include { INTERLEAVE; DEINTERLEAVE; PIGZ_COMPRESS    }  from './modules/base'
include { HUMANN; CHECK_MPA; HUMANN_MERGE }             from './modules/humann'
include { GROOT }                                       from './modules/groot'
include { MULTIQC }                                     from './modules/reports'
/* 
 *   DSL2 allows to reuse channels
 */
reads_ch = Channel
        .fromFilePairs(reads, checkIfExists: true)


 
workflow {

  // Run Kraken2
  KRAKEN2_REPORT(reads_ch, krakenDbPath)
  
  // Check MPA version and run Humann + Metaphlan
  CHECK_MPA(metaphlanDB, params.mpa_ver)
  INTERLEAVE(reads_ch)
  HUMANN( INTERLEAVE.out, chocophlanDB, unirefDB, metaphlanDB, params.mpa_ver, CHECK_MPA.out)

  // TODO: Check Groot DB and read length compatibility
  if (params.grootdb != false) {
    GROOT(INTERLEAVE.out, params.grootdb)
  }
  

  // Make summaries
  HUMANN_MERGE(HUMANN.out.map{it -> it[1]}.collect())
  MULTIQC(KRAKEN2_REPORT.out.report.collect())
}
