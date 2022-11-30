/* 
 
*/

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

include { KRAKEN2_HOST; KRAKEN2_REPORT } from './modules/kraken'
include { INTERLEAVE; DEINTERLEAVE; PIGZ_COMPRESS    }  from './modules/base'
include { HUMANN; CHECK_MPA }                       from './modules/humann'
/* 
 *   DSL2 allows to reuse channels
 */
reads_ch = Channel
        .fromFilePairs(reads, checkIfExists: true)


 
workflow {
  CHECK_MPA(metaphlanDB, params.mpa_ver)

  KRAKEN2_REPORT(reads_ch, krakenDbPath)
  INTERLEAVE(reads_ch)
  HUMANN( INTERLEAVE.out, chocophlanDB, unirefDB, metaphlanDB, params.mpa_ver)
}