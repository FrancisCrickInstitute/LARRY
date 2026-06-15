<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-larry_logo_2_dark.svg">
    <img alt="nf-core/larry" src="docs/images/nf-core-larry_logo_2_light.svg">
  </picture>
</h1>


## Introduction

**nf-core/larry** is a bioinformatics pipeline that determines the clonal composition of cells. Cells have been labeled with (a) LARRY barcode(s). Clones are cells that originate from the same LARRY labelled cel (belong to the same lineage). The pipeline requires 10x Chromium Single Cell data (GEX) and paired-end Illumina sequences of PCR-enriched LARRY labels (LARRY). It takes a samplesheet and FASTQ files as input 

<p align="center">
    <img src="docs/images/workflow_diagram.svg" alt="nf-core/larry schematic workflow" width="70%"
</p>

1. Cellranger count module (GEX) ([`cellranger_count`](https://nf-co.re/modules/cellranger_count/))
2. Obtain the unmapped reads (GEX) ([`samtools_view`](https://nf-co.re/modules/samtools_view/))
3. Obtain reads mapped to the LARRY construct (GEX) ([`samtools_view`](https://nf-co.re/modules/samtools_view/))
4. Merge the unmapped and LARRY reads (GEX) ([`samtools_merge`](https://nf-co.re/modules/samtools_merge/))
5. Sort bam file (GEX) ([`samtools_sort`](https://nf-co.re/modules/samtools_sort/))
6. Convert 10x bam to FASTQ (GEX) ([`bamtofastq10x`](https://nf-co.re/modules/bamtofastq10x/))
7. Concatenate fastq files from same sample (GEX + LARRY) ([`cat_fastq`](https://nf-co.re/modules/cat_fastq/))
8. Remove adapter sequences (GEX + LARRY) ([`cutadapt`](https://nf-co.re/modules/cutadapt/))
9. Check for valid LARRY label (GEX + LARRY) ([`cutadapt`](https://nf-co.re/modules/cutadapt/))
10. Select first 28 nucleotides from the second read to obtain 10x barcode and UMI (CELLUMI) (GEX + LARRY) ([`cutadapt`](https://nf-co.re/modules/cutadapt/))
11. Filter LARRY and CULLUMI reads based on length (GEX + LARRY) ([`cutadapt`](https://nf-co.re/modules/cutadapt/))
12. Obtain cell barcodes from Cellranger that make the treshold (GEX) ([`gunzip`](https://nf-co.re/modules/gunzip/)) ([`cat_cat`](https://nf-co.re/modules/cat_cat/))
13. Extract UMI barcode from a read and add it to the read name (GEX + LARRY) ([`umitools_extract`](https://nf-co.re/modules/umitools_extract/))
14. Run bespoke script to obtain the LARRY label counts, determine cutoff and combine clone ids with multiple LARRY labels. For each sample, outputs the clone assignment table, a QC exclusion table (one row per excluded cell with reason), and a QC summary table (mutually exclusive cell counts per filtering step).

## Usage


First, prepare a samplesheet with your input data that looks as follows:

`input.csv`:

```csv
sample,fastq_1,fastq_2
sample1_GEX_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
sample1_GEX_2,AEG588A2_S1_L002_R1_001.fastq.gz,AEG588A2_S1_L002_R2_001.fastq.gz
sample1_LARRY_1,AEG588A3_S1_L002_R1_001.fastq.gz,AEG588A3_S1_L002_R2_001.fastq.gz
```

Each row represents a pair of fastq files (paired end).
Be aware! For the LARRY libraries: The first read needs to be the LARRY barcode and the second read the CELL + UMI barcode.

The sample name should be structured as follows:
1. Sample name
2. GEX/LARRY (depending on the library)
3. File number. Files with same sample name and GEX/LARRY should be enumerated: 1, 2, 3 ...

These three elements should be connected by an "_"

Clone the nextflow pipeline repository, go into the directory and switch to the dev branch.

```bash
git clone git@github.com:FrancisCrickInstitute/LARRY.git
cd LARRY/
git switch dev
```

The pipeline runs with Singularity. The Singularity module needs to be available.

Being in the nextflow pipeline directory, run the pipeline like this:

```bash
nextflow main.nf \
   --input <SAMPLESHEET> \
   --outdir <OUTDIR> \
   --cellranger_reference <CELLRANGER REFERENCE (absolute path)> \
   -w <SCRATCH WORK DIR> \
   -resume
```

`-w` sets the directory where Nextflow writes intermediate files (recommended: a scratch area). `-resume` allows restarting from the last successful step if the run fails.

## Parameters

These parameters can be set in the nextflow.config file

- ham_larry: hamming distance between LARRY sequences to be considered the same (default 3)
- ham_umi: hamming distance between UMI's to be considered the same (default 1)
- gex_and_larry: true if a GEX and LARRY library is present for the samples, false if only a GEX library is present (default true)
- epsilon: Parameter to determine the LARRY expression level cutoff. The lowest counts are stepwise removed and shifted by -1 until relative change in count mean =< epsilon (default 0.01)
- jaccard_cutoff: clones share many of the same cells, jaccard distance between clones >= jaccard_clone, are merged (default 0.5)

## Pipeline output

For every sample there is an output for the GEX and LARRY library, if both were available.
Every library has four files:

1. `*_clone_output.csv` — clone assignment table, 3 columns:
    - Clone: Name of the clone that the cell belongs to.
    - Cell: Name of the cell.
    - Count: Number of LARRY reads counted within the cell.
2. `*.png` — figure depicting the LARRY signal cutoff.
    - X axis: number of reads for a particular LARRY label within a particular cell.
    - Y axis: frequency of occurrence of that particular read count (log scale).
    - Red line: the mean count after stepwise removal and shift of the lowest counts. When the slope flattens the cutoff is set (controlled by the `epsilon` parameter).
    - Blue vertical line: signal cutoff. LARRY signals below this count are not considered.
3. `*_qc_exclusions.csv` — one row per excluded cell barcode, with columns `sample_id`, `step`, and `cell_barcode`. Steps are: `ambiguous_LARRY_barcode`, `low_UMI_count`, `ambiguous_clone_assignment`.
4. `*_qc_summary.csv` — one row per filtering step, with columns `sample_id`, `step`, `n_umi_groups_removed`, and `n_excluded_cell_barcodes`. Counts are mutually exclusive across steps and restricted to CellRanger-called cells. A `total_excluded` row is included.

## Credits

This was originally written by Jasper Depotter.

We thank the following people for their extensive assistance in the development of this pipeline:

- Stephanie Strohbuecker
- Giulia Boezio

## Test

This test data is not in the github repository yet!

To test if the complete pipeline runs succesfully.
Run in the nextflow pipeline directory:

```bash
nf-test test tests/main.nf.test
```
