nextflow.enable.moduleBinaries = true

include { FETCH_READS }                        from './modules/fetch_reads.nf'
include { TAG_BARCODES }                       from './modules/tag_barcodes/main.nf'
include { PREPARE_REFERENCE ; BWAMEM2_INDEX }  from './modules/reference.nf'
include { ALIGN ; MERGE_MARKDUP }              from './modules/align.nf'
include { BARCODE_MATRIX }                     from './modules/barcode_matrix/main.nf'
include { MAKE_BCM }                           from './modules/make_bcm/main.nf'

params {
    // CSV with a header and the columns sample,accession (SRA run accession).
    input: Path

    // Reference FASTA: plain, gzip or bgzip.
    reference: Path

    // Chromosomes to build tables for: reference sequence names matching this regex.
    chromosomes: String = '^Cec[0-9]+$'

    // Read pairs to take from the start of each run; 0 downloads the full run.
    max_pairs: Integer = 0

    // Fraction of barcodes kept (all read pairs of a kept barcode are kept).
    barcode_fraction: Float = 0.25
    seed: Integer = 1

    // Read pairs per alignment chunk.
    chunk_pairs: Integer = 10000000

    // Bin size in bp of the co-occurrence tables.
    bin_size: Integer = 5000

    // Also publish the duplicate-removed BAMs.
    publish_bam: Boolean = false
}

workflow {
    main:
    samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.accession) }

    reads  = FETCH_READS(samples, params.max_pairs)

    ref   = PREPARE_REFERENCE(params.reference)
    index = BWAMEM2_INDEX(ref.fasta)
}
