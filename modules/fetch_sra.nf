// Alternative to FETCH_READS: take the run from NCBI instead of ENA.
//
// The .sra object is downloaded from the location NCBI reports for the accession (a public
// cloud bucket), checked against its md5 sum, and converted to paired FASTQ with one
// fastq-dump per CPU, each on its own range of spots. The parts are concatenated in order
// (gzip members can be concatenated). max_pairs > 0 converts only the first max_pairs spots.
process FETCH_SRA {
    tag "${sample}"
    container 'quay.io/biocontainers/sra-tools:3.2.1--h4304569_1'

    input:
    tuple val(sample), val(accession)
    val max_pairs

    output:
    tuple val(sample), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    """
    set -euo pipefail
    loc=\$(curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=${accession}")
    url=\$(echo "\$loc" | grep -o '"link": *"[^"]*"' | head -n 1 | cut -d'"' -f4)
    md5=\$(echo "\$loc" | grep -o '"md5": *"[^"]*"' | head -n 1 | cut -d'"' -f4)
    [ -n "\$url" ] && [ -n "\$md5" ] || { echo "ERROR: NCBI reports no download location for ${accession}: \$loc" >&2; exit 1; }

    t0=\$(date +%s)
    curl -fsSL --retry 10 --retry-delay 15 --retry-all-errors -o ${accession}.sra "\$url"
    t1=\$(date +%s)
    bytes=\$(stat -c %s ${accession}.sra)
    awk -v n=\$bytes -v s=\$((t1 - t0)) 'BEGIN{printf "[FETCH_SRA] download: %.1f GB in %d s, %.1f MB/s\\n", n / 1e9, s, n / 1e6 / (s > 0 ? s : 1)}' >&2
    echo "\$md5  ${accession}.sra" | md5sum -c - >&2
    t2=\$(date +%s)

    spots=\$(vdb-dump --id_range ${accession}.sra | sed 's/.*row-count = //' | tr -cd '0-9')
    n=\$spots
    if [ ${max_pairs} -gt 0 ] && [ ${max_pairs} -lt \$spots ]; then n=${max_pairs}; fi
    parts=${task.cpus}
    size=\$(( (n + parts - 1) / parts ))

    pids=""
    for i in \$(seq 0 \$((parts - 1))); do
        first=\$((i * size + 1))
        last=\$(( (i + 1) * size ))
        [ \$last -gt \$n ] && last=\$n
        [ \$first -gt \$last ] && continue
        part=\$(printf 'part_%03d' \$i)
        fastq-dump -N \$first -X \$last --split-files --gzip \\
            --defline-seq '@\$ac.\$si/\$ri' --defline-qual '+' \\
            -O \$part ${accession}.sra > \$part.log 2>&1 &
        pids="\$pids \$!"
    done
    for pid in \$pids; do
        wait \$pid || { cat part_*.log >&2; echo "ERROR: fastq-dump failed for ${accession}" >&2; exit 1; }
    done

    written=\$(cat part_*.log | awk '\$1 == "Written" {s += \$2} END {print s + 0}')
    [ "\$written" -eq "\$n" ] || { echo "ERROR: wrote \$written of \$n spots for ${accession}" >&2; exit 1; }
    cat part_*/${accession}_1.fastq.gz > ${accession}_1.fastq.gz
    cat part_*/${accession}_2.fastq.gz > ${accession}_2.fastq.gz
    t3=\$(date +%s)
    awk -v n=\$n -v s=\$((t3 - t2)) -v p=\$parts 'BEGIN{printf "[FETCH_SRA] convert: %d spots in %d s with %d processes, %.0f spots/s\\n", n, s, p, n / (s > 0 ? s : 1)}' >&2
    rm -rf part_* ${accession}.sra
    ls -l ${accession}_1.fastq.gz ${accession}_2.fastq.gz >&2
    """

    stub:
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${accession}_1.fastq.gz
    printf '@r/2\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${accession}_2.fastq.gz
    """
}
