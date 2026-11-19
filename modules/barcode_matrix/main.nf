// Count, for every pair of bins on one chromosome, the barcodes seen in both bins.
// Output: binA_start, binB_start, count; symmetric, diagonal once.
process BARCODE_MATRIX {
    tag "${sample}.${chrom}"
    conda 'bioconda::samtools=1.21 conda-forge::python=3.12 conda-forge::numpy=2.1 conda-forge::scipy=1.14 conda-forge::pandas=2.2'

    input:
    tuple val(sample), path(bam), path(bai), val(chrom)
    val bin_size
    val bin_label

    output:
    tuple val(sample), val(chrom), path("${sample}.${chrom}_${bin_label}_bin.stat.gz"), emit: stat
    tuple val(sample), val(chrom), path("${sample}.${chrom}_${bin_label}_bin.summary.json"), emit: summary

    script:
    """
    barcode_matrix.sh ${bam} ${chrom} ${bin_size} \\
        ${sample}.${chrom}_${bin_label}_bin.stat.gz \\
        ${sample}.${chrom}_${bin_label}_bin.summary.json \\
        ${task.cpus}
    """

    stub:
    """
    printf '0\\t0\\t5\\n' | gzip -c > ${sample}.${chrom}_${bin_label}_bin.stat.gz
    echo '{"stub": true}' > ${sample}.${chrom}_${bin_label}_bin.summary.json
    """
}
