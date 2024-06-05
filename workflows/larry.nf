/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                             } from '../modules/nf-core/fastqc/main'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main'
include { CUTADAPT as CUTADAPT_remove_adapt  } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_valid_larry   } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_cut_umi       } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_len_filter    } from '../modules/nf-core/cutadapt/main'
include { UMITOOLS_EXTRACT } from '../modules/nf-core/umitools/extract/main'
include { paramsSummaryMap                   } from 'plugin/nf-validation'
include { paramsSummaryMultiqc               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_larry_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LARRY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



    //
    //Module concatenate
    //

    /*
    CONCAT_FASTQ(
        CONCAT_FASTQ
    )
    */

    //
    //Module cutadapt
    //

    CUTADAPT_remove_adapt(
        ch_samplesheet
    )

    

    CUTADAPT_valid_larry(
        CUTADAPT_remove_adapt.out.reads
    )

    

    Channel
        CUTADAPT_valid_larry.out.reads
                            .map{ reads ->
                            def meta = reads[0]
                            meta.single_end = true

                            def secondRead = reads[1][1]

                            return new Tuple(meta,secondRead)

                            }
                            .set{CUTADAPT_cut_umi_input}


    Channel
        CUTADAPT_valid_larry.out.reads
                            .map{ reads ->
                            def meta = reads[0]

                            def firstRead = reads[1][0]

                            [["id" : meta.id , "read" : firstRead]]
                            }
                            .collect()
                            .collectEntries(){item ->
                                            [(item.id): item.read]
                            }
                            .set{firstRead}
    
    CUTADAPT_cut_umi(
        CUTADAPT_cut_umi_input
    )

    Channel
        CUTADAPT_cut_umi.out.reads
                        .map{ reads ->
                        def meta = reads[0]
                        meta.single_end = false

                        def shortRead = reads[1]

                        def first_read_el = firstRead[meta.id].value

                        updatedReads = [first_read_el , shortRead]

                        return new Tuple(meta,updatedReads)

                        }
                        .set{CUTADAPT_len_filter_input}

    CUTADAPT_len_filter(
        CUTADAPT_len_filter_input
    )
                        

    /*
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        CUTADAPT_cut_umi.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    */

    emit:
    cutadapt_ra = CUTADAPT_len_filter.out.reads
    //multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
