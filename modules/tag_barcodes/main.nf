// Move the stLFR barcode from the tail of read 2 into the read name (name#A_B_C), trim
// read 2 to the 100 bp insert, drop pairs without a valid barcode, optionally keep a
// random fraction of barcodes, and write the result as FASTQ chunks for parallel alignment.
// The script is a module binary (resources/usr/bin).
process TAG_BARCODES {
    tag "${sample}"
    conda 'conda-forge::python=3.12 conda-forge::pigz=2.8'

    input:
    tuple val(sample), path(r1), path(r2)
    path whitelist
    val chunk_pairs
    val barcode_fraction
    val seed

    output:
    tuple val(sample), path("${sample}.chunk*_{1,2}.fastq.gz"), emit: chunks
    tuple val(sample), path("${sample}.tag_stats.json"), emit: stats

    script:
    """
    stlfr_tag_barcodes.py \\
        --r1 ${r1} \\
        --r2 ${r2} \\
        --whitelist ${whitelist} \\
        --out-prefix ${sample} \\
        --chunk-pairs ${chunk_pairs} \\
        --barcode-fraction ${barcode_fraction} \\
        --seed ${seed} \\
        --stats ${sample}.tag_stats.json \\
        --threads ${task.cpus}
    """

    stub:
    """
    for c in 0000 0001; do
        printf '@r#1_2_3/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${sample}.chunk\${c}_1.fastq.gz
        printf '@r#1_2_3/2\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${sample}.chunk\${c}_2.fastq.gz
    done
    echo '{"stub": true}' > ${sample}.tag_stats.json
    """
}
