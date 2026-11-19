// Decompress the reference (plain, gzip or bgzip) and index it with samtools faidx.
process PREPARE_REFERENCE {
    conda 'bioconda::samtools=1.21'

    input:
    path fasta

    output:
    path 'genome.fa', emit: fasta
    path 'genome.fa.fai', emit: fai

    script:
    """
    set -euo pipefail
    if gzip -t ${fasta} 2>/dev/null; then
        gzip -dc ${fasta} > genome.fa
    else
        cp ${fasta} genome.fa
    fi
    samtools faidx genome.fa
    """

    stub:
    """
    touch genome.fa
    printf 'Cec1\\t100000\\t6\\t80\\t81\\nCec2\\t80000\\t6\\t80\\t81\\nCec_scaffold1\\t5000\\t6\\t80\\t81\\n' > genome.fa.fai
    """
}

process BWAMEM2_INDEX {
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21'

    input:
    path fasta

    output:
    path 'bwamem2'

    script:
    """
    mkdir bwamem2
    bwa-mem2 index -p bwamem2/genome ${fasta}
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/genome.0123 bwamem2/genome.amb bwamem2/genome.ann bwamem2/genome.bwt.2bit.64 bwamem2/genome.pac
    """
}
