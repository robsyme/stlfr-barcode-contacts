// Alignment statistics of the duplicate-removed BAM, in the formats MultiQC reads.
process BAM_STATS {
    tag "${sample}"
    conda 'bioconda::samtools=1.21'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.flagstat"), path("${sample}.stats")

    script:
    """
    set -euo pipefail
    samtools flagstat -@ ${task.cpus} ${bam} > ${sample}.flagstat
    samtools stats -@ ${task.cpus} ${bam} > ${sample}.stats
    """

    stub:
    """
    touch ${sample}.flagstat ${sample}.stats
    """
}
