/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CELLRANGER_COUNT } from '../modules/nf-core/cellranger/count/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_unmapped} from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_larry} from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { BAMTOFASTQ10X } from '../modules/nf-core/bamtofastq10x/main'
include { CAT_FASTQ} from '../modules/nf-core/cat/fastq/main'
include { CUTADAPT as CUTADAPT_remove_adapt  } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_valid_larry   } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_cut_umi       } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_len_filter    } from '../modules/nf-core/cutadapt/main'
include { GUNZIP                   } from '../modules/nf-core/gunzip/main'
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
include { UNIQUE } from '../modules/local/unique/main'
include { UMITOOLS_EXTRACT                   } from '../modules/nf-core/umitools/extract/main'
include { FASTQ_SIZE } from '../modules/local/fastq_size/main'
include { GATHER_BARCODE                   } from '../modules/local/gather_barcode/main'
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
    //Filter 10x data
    //


    ch_samplesheet
        .filter { it[0].id.contains('GEX') }
        .set{ch_samplesheet_gex}


    //
    //Module cellranger count
    //

    CELLRANGER_COUNT(ch_samplesheet_gex , params.cellranger_reference)

    //
    //Extract BAM files
    //

    CELLRANGER_COUNT.out.outs
                    .map{meta , fles -> tuple(meta , fles[-7] , fles[-6])}
                    .set{bam_files}

    //
    //Module samtool views for unmapped reads
    //

    SAMTOOLS_VIEW_unmapped(bam_files , tuple([:] , []) , [])

    
    //
    //Module samtool views for LARRY reads
    //

    SAMTOOLS_VIEW_larry(bam_files , tuple([:] , []) , [])

    //
    //Prepare bamfiles to merge
    //

    SAMTOOLS_VIEW_unmapped.out.bam
                        .concat(SAMTOOLS_VIEW_larry.out.bam)
                        .map{ meta, fastq -> 
                            def lastUnderscoreIndex = meta.id.lastIndexOf('_')
                            def modifiedMetaId = meta.id.substring(0, lastUnderscoreIndex)
                            tuple( modifiedMetaId ,  fastq )}
                        .groupTuple()
                        .map{ meta, fastq -> tuple( [id : meta , single_end : false], fastq.flatten() )}
                        .set{bam_file_tomerge}

    //
    //Module samtool merge
    //

    SAMTOOLS_MERGE(bam_file_tomerge , tuple([:] , []) , tuple([:] , []))

    //
    //Module samtool sort
    //

    SAMTOOLS_SORT(SAMTOOLS_MERGE.out.bam , tuple([:] , []))

    //
    //Module BAMTOFASTQ10X 
    //

    BAMTOFASTQ10X(SAMTOOLS_SORT.out.bam)


    //
    //Combine LARRY and 10x data again
    //

    ch_samplesheet.filter { it[0].id.contains('LARRY')}
                    .map{ meta, fastq -> 
                        def lastUnderscoreIndex = meta.id.lastIndexOf('_')
                        def modifiedMetaId = meta.id.substring(0, lastUnderscoreIndex)
                        tuple( modifiedMetaId ,  fastq )}
                    .groupTuple()
                    .map{ meta, fastq -> tuple( [id : meta , single_end : false], fastq.flatten() )}
                    .concat(BAMTOFASTQ10X.out.fastq)
                    .set{gex_larry_together}

    //
    //Module concatenate
    //

    CAT_FASTQ(
        gex_larry_together
    )

    
    //
    //Reverse R1 and R2 order
    //


    CAT_FASTQ.out.reads
                .map{ meta, fastq -> tuple(meta , [fastq[1] , fastq[0]])}
                .set{gex_larry_together_reverse}


    //
    //Module cutadapt remove adapter
    //

    CUTADAPT_remove_adapt(
        gex_larry_together_reverse
    )

    
    //
    //Module cutadapt check valid larry
    //

    CUTADAPT_valid_larry(
        CUTADAPT_remove_adapt.out.reads
    )

    
    //
    //Convert to a single-end data with only the cell and umi barcode.
    //

    CUTADAPT_valid_larry.out.reads
                        .map{ meta, fastq -> 
                            meta.single_end  = true
                            return tuple(meta , fastq[1])}
                        .set{CUTADAPT_cut_umi_input}

    //
    //Link sample identity to the LARRY barcode
    //

    CUTADAPT_valid_larry.out.reads
                        .map{ meta, fastq -> 
                        meta.single_end  = true
                        return tuple(meta, fastq[0])}
                        .set{firstRead}
    
    //
    //Module cutadapt select first 28 nucleotides of cellumi reads
    //

    CUTADAPT_cut_umi(
        CUTADAPT_cut_umi_input
    )

    //
    //Combine LARRY and CELLUMI reads again.
    //

    firstRead.concat(CUTADAPT_cut_umi.out.reads)
                .map{ meta, fastq -> tuple( meta.id, fastq )}
                .groupTuple()
                .map{ meta, fastq -> tuple( [id : meta , single_end : false], fastq.flatten() )}
                .set{CUTADAPT_len_filter_input}


    //
    //Module cutadapt filter LARRY and CELLUMI reads based on length.
    //

    CUTADAPT_len_filter(
        CUTADAPT_len_filter_input
    )


    //
    //LARRY and CELLUMI need to change from position if we want to use UMITOOLS_EXTRACT
    //

    CUTADAPT_len_filter.out.reads
                        .map{meta , fastq -> tuple(meta.id , [fastq[1],fastq[0]])}
                        .set{CUTADAPT_output}

    //
    //Get barcode files
    //

    CELLRANGER_COUNT.out.outs
                    .map{meta , fles -> tuple(meta, fles[-4])}
                    .set{barcodes}

    //
    //Gunzip the barcodes
    //

    GUNZIP(barcodes)


    //
    //Prepare files to concatenate
    //


    GUNZIP.out.gunzip
                .map{ meta, barcode -> 
                    def lastUnderscoreIndex = meta.id.lastIndexOf('_')
                    def modifiedMetaId = meta.id.substring(0, lastUnderscoreIndex)
                    tuple( modifiedMetaId ,  barcode )}
                .groupTuple()
                .map{ meta, barcode -> tuple( [id : meta , single_end : false], barcode.flatten() )}
                .set{barcodes_cat}

    
    //
    //Combine barcode files frome same sample
    //

    CAT_CAT(barcodes_cat)


    //
    //Make barcodes unique and remove "-1" suffix
    //

    UNIQUE(CAT_CAT.out.file_out)

    //
    //Add whitelist to read files
    //

    UNIQUE.out.file_out
            .map{meta , file_path -> tuple(meta.id, file_path)}
            .set{barcodes_gex}

    if (params.gex_and_larry){

        UNIQUE.out.file_out
        .map{meta , file_path -> tuple(meta.id.substring(0, meta.id.lastIndexOf('_')) + "_LARRY", file_path)}
        .concat(barcodes_gex)
        .set{barcodes_gex_larry}

    }

    else {
        barcodes_gex.set{barcodes_gex_larry}

    }

    CUTADAPT_output.concat(barcodes_gex_larry)
                            .groupTuple()
                            .map{meta , elmnts -> tuple([id : meta , single_end : false], elmnts[0] , elmnts[1] )}
                            .set{UMITOOLS_EXTRACT_input}

    //
    //Run UMITOOLS script
    //

    UMITOOLS_EXTRACT(
        UMITOOLS_EXTRACT_input
    )


    UMITOOLS_EXTRACT.out.reads
        .map{ meta, fastq -> tuple( [id : meta.id , single_end : true] , fastq[1] )}
        .set{COUNT_FASTQ_input}


    FASTQ_SIZE(COUNT_FASTQ_input)

    FASTQ_SIZE.out.outs
        .filter{v -> v[2].toInteger() > 1}
        .map{meta , fastq , file_size -> tuple(meta , fastq)}
        .set{GATHER_BARCODE_input}

    //
    //Gather the barcode: at some point I need to implement that LARRY and 10X library are combined together
    //
    
    GATHER_BARCODE(
        GATHER_BARCODE_input
        )


    emit:
    larry_barcodes = GATHER_BARCODE.out.outs
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
