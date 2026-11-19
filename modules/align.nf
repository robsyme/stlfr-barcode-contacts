// Align one FASTQ chunk. bwa-mem2 drops the trailing /1 and /2, so the BAM read names
// are "id#A_B_C". fixmate -m adds the mate score tags that markdup needs later.
process ALIGN {
    tag "${sample}.${chunk}"
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21'

    input:
    tuple val(sample), val(chunk), path(r1), path(r2)
    path index

    output:
    tuple val(sample), path("${sample}.${chunk}.sorted.bam")

    script:
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    """
    set -euo pipefail
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample}.${chunk}\\tSM:${sample}\\tLB:${sample}\\tPL:DNBSEQ' \\
        ${index}/genome ${r1} ${r2} \\
        | samtools fixmate -m -u - - \\
        | samtools sort -@ ${sort_threads} -m 1G -T sort_tmp -o ${sample}.${chunk}.sorted.bam -
    """

    stub:
    """
    touch ${sample}.${chunk}.sorted.bam
    """
}

// Merge the chunk BAMs of a sample and remove duplicates.
process MERGE_MARKDUP {
    tag "${sample}"
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21'

    input:
    tuple val(sample), path(bams)

    output:
    tuple val(sample), path("${sample}.markdup.bam"), path("${sample}.markdup.bam.bai"), emit: bam
    tuple val(sample), path("${sample}.markdup_stats.txt"), emit: stats

    script:
    """
    set -euo pipefail
    samtools merge -@ ${task.cpus} -u -o - ${bams} \\
        | samtools markdup -r -@ ${task.cpus} -T markdup_tmp -s -f ${sample}.markdup_stats.txt - ${sample}.markdup.bam
    samtools index -@ ${task.cpus} ${sample}.markdup.bam
    """

    stub:
    """
    touch ${sample}.markdup.bam ${sample}.markdup.bam.bai ${sample}.markdup_stats.txt
    """
}
