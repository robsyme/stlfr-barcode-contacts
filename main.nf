nextflow.enable.moduleBinaries = true

include { FETCH_READS }                        from './modules/fetch_reads.nf'
include { FETCH_SRA }                           from './modules/fetch_sra.nf'
include { TAG_BARCODES }                       from './modules/tag_barcodes/main.nf'
include { PREPARE_REFERENCE ; BWAMEM2_INDEX }  from './modules/reference.nf'
include { ALIGN ; MERGE_MARKDUP }              from './modules/align.nf'
include { BARCODE_MATRIX }                     from './modules/barcode_matrix/main.nf'
include { MAKE_BCM }                           from './modules/make_bcm/main.nf'
include { BAM_STATS }                          from './modules/bam_stats.nf'
include { MULTIQC }                            from './modules/multiqc/main.nf'

params {
    // CSV with a header and the columns sample,accession (SRA run accession).
    input: Path

    // Reference FASTA: plain, gzip or bgzip.
    reference: Path

    // Chromosomes to build tables for: reference sequence names matching this regex.
    chromosomes: String = '^Cec[0-9]+$'

    // Where to download reads from: 'ena' (FASTQ files) or 'sra' (NCBI .sra object, converted).
    read_source: String = 'ena'

    // Read pairs to take from the start of each run; 0 downloads the full run.
    max_pairs: Integer = 0

    // Fraction of barcodes kept (all read pairs of a kept barcode are kept).
    barcode_fraction: Float = 0.25
    seed: Integer = 1

    // Read pairs per alignment chunk.
    chunk_pairs: Integer = 10000000

    // Bin size in bp of the co-occurrence tables.
    bin_size: Integer = 5000

    // stLFR barcode whitelist; defaults to the copy shipped with the pipeline.
    barcode_whitelist: Path = "${projectDir}/modules/tag_barcodes/resources/barcode_whitelist.tsv"

    // Also publish the duplicate-removed BAMs.
    publish_bam: Boolean = false
}

workflow {
    main:
    samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.accession) }

    reads  = params.read_source == 'sra'
        ? FETCH_SRA(samples, params.max_pairs)
        : FETCH_READS(samples, params.max_pairs)
    tagged = TAG_BARCODES(reads, params.barcode_whitelist, params.chunk_pairs, params.barcode_fraction, params.seed)

    // One element per chunk: sample, chunk id, read 1, read 2.
    chunks = tagged.chunks.flatMap { sample, files ->
        def list = files instanceof List ? files : [files]
        list
            .groupBy { f -> (f.name =~ /\.chunk(\d+)_[12]\.fastq\.gz$/)[0][1] }
            .collect { id, pair ->
                def sorted = pair.sort { f -> f.name }
                tuple(sample, "chunk${id}", sorted[0], sorted[1])
            }
    }

    ref   = PREPARE_REFERENCE(params.reference)
    index = BWAMEM2_INDEX(ref.fasta)

    aligned = ALIGN(chunks, index)

    // Group the chunk BAMs by sample. The group size is known from the tagging output, so a
    // sample is merged as soon as its own chunks are aligned. BAMs are sorted by name to keep
    // the task hash stable across resumed runs.
    n_chunks = tagged.chunks.map { sample, files -> tuple(sample, (files instanceof List ? files.size() : 1).intdiv(2)) }
    grouped  = aligned
        .combine(n_chunks, by: 0)
        .map { sample, bam, n -> tuple(groupKey(sample, n), bam) }
        .groupTuple()
        .map { key, bams -> tuple(key.target, bams.sort { b -> b.name }) }
    merged   = MERGE_MARKDUP(grouped)

    chroms = ref.fai
        .splitCsv(sep: '\t')
        .map { row -> row[0] }
        .filter { name -> name ==~ params.chromosomes }

    bin_label = params.bin_size % 1000 == 0 ? "${params.bin_size.intdiv(1000)}kb" : "${params.bin_size}bp"
    matrix  = BARCODE_MATRIX(merged.bam.combine(chroms), params.bin_size, bin_label)
    bundles = MAKE_BCM(matrix.stat, ref.fai, params.bin_size)

    // QC report: tagging, duplicate, alignment and table statistics of all samples.
    // Files are sorted by name so the report task has a stable hash across resumed runs.
    bam_stats = BAM_STATS(merged.bam)
    qc_files  = tagged.stats.map { _sample, f -> f }
        .mix(merged.stats.map { _sample, f -> f })
        .mix(bam_stats.flatMap { _sample, flagstat, stats -> [flagstat, stats] })
        .mix(matrix.summary.map { _sample, _chrom, f -> f })
        .toSortedList { a, b -> a.name <=> b.name }
    report    = MULTIQC(qc_files, file("${projectDir}/assets/multiqc_config.yaml"))

    publish:
    multiqc_report = report.report
    multiqc_data   = report.data
    bam_stats      = bam_stats.map { _sample, flagstat, stats -> [flagstat, stats] }
    bcm            = bundles.map { sample, _chrom, f -> tuple(sample, f) }
    stat           = matrix.stat.map { sample, _chrom, f -> tuple(sample, f) }
    matrix_summary = matrix.summary.map { sample, _chrom, f -> tuple(sample, f) }
    tag_stats      = tagged.stats.map { _sample, f -> f }
    markdup_stats  = merged.stats.map { _sample, f -> f }
    bam            = params.publish_bam ? merged.bam : channel.empty()
}

output {
    bcm {
        path { sample, _f -> "bcm/${sample}" }
    }
    stat {
        path { sample, _f -> "stat/${sample}" }
    }
    matrix_summary {
        path { sample, _f -> "stat/${sample}" }
    }
    tag_stats {
        path 'tag_stats'
    }
    markdup_stats {
        path 'markdup_stats'
    }
    multiqc_report {
        path 'multiqc'
    }
    multiqc_data {
        path 'multiqc'
    }
    bam_stats {
        path 'bam_stats'
    }
    bam {
        path { sample, _bam, _bai -> "bam/${sample}" }
    }
}
