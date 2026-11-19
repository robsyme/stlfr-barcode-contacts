// Build the viewer bundle for one chromosome: contacts plus the hotspot and junction calls.
process MAKE_BCM {
    tag "${sample}.${chrom}"
    conda 'conda-forge::python=3.12 conda-forge::numpy=2.1 conda-forge::scipy=1.14'

    input:
    tuple val(sample), val(chrom), path(stat)
    path fai
    val bin_size

    output:
    tuple val(sample), val(chrom), path("${chrom}.bcm")

    script:
    """
    preprocess.py --input ${stat} --name ${chrom} --fai ${fai} --bin ${bin_size} --out .
    """

    stub:
    """
    touch ${chrom}.bcm
    """
}
