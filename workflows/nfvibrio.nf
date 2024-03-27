/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNfvibrio.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_TRIM_FASTP_FASTQC } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTP                     } from '../modules/nf-core/fastp/main'
include { CAT_FASTQ                     } from '../modules/nf-core/cat/fastq/main'
include { TRIMMOMATIC                   } from '../modules/nf-core/trimmomatic/main'
include { CONFINDR                      } from '../modules/local/confindr/main'
include { KRAKEN2_KRAKEN2               } from '../modules/nf-core/kraken2/kraken2/main'
include { SHOVILL                       } from '../modules/nf-core/shovill/main'
include { SNIPPY_RUN                    } from '../modules/nf-core/snippy/run/main'
include { ABRICATE_RUN                  } from '../modules/nf-core/abricate/run/main'
include { ROARY                         } from '../modules/nf-core/roary/main'
include { RENAME                         } from '../modules/local/rename.nf'
//include { GUBBINS                       } from '../modules/nf-core/gubbins/main'
include { QUAST                         } from '../modules/nf-core/quast/main'
include { PROKKA                        } from '../modules/nf-core/prokka/main'
include { AMRFINDERPLUS_UPDATE          } from '../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN             } from '../modules/nf-core/amrfinderplus/run/main'
include { PLOT_VF_HEATMAP               } from '../modules/local/plot/vf_heatmap/main'
include { MLST                          } from '../modules/nf-core/mlst/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NFVIBRIO {

    ch_versions = Channel.empty()
    
    ch_reference_fasta          = Channel.fromPath(params.reference_fasta).first()
    ch_reference_gff            = Channel.fromPath(params.reference_gff).first()
    ch_global_meta              = Channel.from([ id: "main", single_end:false ]).first()
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SUBWORKFLOW: Read in samplesheet, validate and stage input files
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Concatenate FastQ files from same sample if required
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Trimmomatic
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    TRIMMOMATIC (
        ch_cat_fastq
    )
    
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SUBWORKFLOW: Read QC and trim adapters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    //ch_trimmed_reads.subscribe{println "value: ${it[0].id}"}
    FASTQ_TRIM_FASTP_FASTQC (
        ch_trimmed_reads,
        Channel.fromPath(params.adapter_path).first(),
        params.save_trimmed_fail,
        false, 
        false,
        false,
      )
    ch_filter_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Kraken2 metagenomic classifier 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_kraken2_multiqc = Channel.empty()
    
    if (params.run_kraken2) {
        ch_kraken2_database = Channel.fromPath(params.kraken2_database)
        KRAKEN2_KRAKEN2 (
            ch_filter_fastq,
            ch_kraken2_database,
            false,  
            true
    )
    ch_kraken2_multiqc = KRAKEN2_KRAKEN2.out.report.map{ it -> it[1]}
    ch_versions        = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first().ifEmpty(null))
    }
    

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Run contamination detection
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    if (params.run_confindr) {
        ch_confindr_fastq  = ch_filter_fastq

        CONFINDR (
            ch_confindr_fastq.collect{it[1]}
        )
        
        ch_versions            = ch_versions.mix(CONFINDR.out.versions)

    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Assembling reads
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    SHOVILL (
        ch_filter_fastq
    )

    ch_versions             = ch_versions.mix(SHOVILL.out.versions)

    RENAME (
        SHOVILL.out.contigs,
        "_contigs.fa"
    )

    ch_contigs              = RENAME.out
    


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Abricate for virulence factors
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ABRICATE_RUN(
        ch_contigs
    )
    ch_versions               = ch_versions.mix(ABRICATE_RUN.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Custom plotting module to generate HTML heatmap
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_abricate_vfs           = ABRICATE_RUN.out.report.collect{it[1]}
    PLOT_VF_HEATMAP(
        ch_abricate_vfs
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Assembly QC 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_quast_contigs        = ch_global_meta.combine(ch_contigs.collect{it[1]}).map{it -> [it[0], it[1..-1]]}
    ch_quast_fasta          = ch_global_meta.combine(ch_reference_fasta)
    ch_quast_gff            = ch_global_meta.combine(ch_reference_gff)

    QUAST(
        ch_quast_contigs,
        ch_quast_fasta,
        ch_quast_gff,
    )
    ch_versions             = ch_versions.mix(QUAST.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Snippy for variant calling & consensus generation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    SNIPPY_RUN(
        ch_filter_fastq,
        ch_reference_fasta
    )
    ch_versions     = ch_versions.mix(SNIPPY_RUN.out.versions)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Prokka genomic annotation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    PROKKA(
        ch_contigs,
        [],
        []
    )
    ch_versions             = ch_versions.mix(PROKKA.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Roary pan-genome generation 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_roary_gffs         = ch_global_meta.combine(PROKKA.out.gff.collect{it[1]}).map{it -> [it[0], it[1..-1]]}
    ROARY(
        ch_roary_gffs
    )
    ch_versions             = ch_versions.mix(ROARY.out.versions)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Gubbins recombination detection
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // ch_gubbins_contigs         = SHOVILL.out.contigs.collect{it[1]}
    // GUBBINS(
    //     ch_gubbins_contigs
    // )
    // ch_versions             = ch_versions.mix(GUBBINS.out.versions)



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULES: Antibiotic resistance
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    AMRFINDERPLUS_UPDATE(
        
    )
    ch_amr_db               = AMRFINDERPLUS_UPDATE.out.db
    ch_versions             = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)

    
    AMRFINDERPLUS_RUN(
        ch_contigs,
        ch_amr_db
    )
    
    ch_versions             = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)
    ch_versions             = ch_versions.mix(AMRFINDERPLUS_RUN.out.tool_version)
    ch_versions             = ch_versions.mix(AMRFINDERPLUS_RUN.out.db_version)

   
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULES: MLST
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    MLST(
        ch_contigs
    )
    ch_mlst_concat          = MLST.out.tsv
    ch_versions             = ch_versions.mix(MLST.out.versions)
   

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: Software and database versions
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
       ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE: MultiQC for generating quality reports
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    
    workflow_summary    = WorkflowNfvibrio.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNfvibrio.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.map{ it -> it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.trim_json.map{ it -> it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(ch_kraken2_multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(SNIPPY_RUN.out.txt.map{ it -> it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.map{ it -> it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
