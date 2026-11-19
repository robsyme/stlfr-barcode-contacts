// Download the paired FASTQ files of one SRA run from ENA.
//
// max_pairs > 0: stream only the first max_pairs read pairs of each file (cheap test runs).
// max_pairs = 0: download both files in full and check them against the ENA md5 sums.
//
// Pipe exit status is not trusted (head closing the pipe kills the upstream commands), so
// the subset mode validates the line count of what was written instead.
process FETCH_READS {
    tag "${sample}"
    // Prebuilt image that carries curl, gzip, awk, coreutils and fastq-dl.
    container 'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3'

    input:
    tuple val(sample), val(accession)
    val max_pairs

    output:
    tuple val(sample), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    def n_lines = (max_pairs as long) * 4
    """
    set -uo pipefail
    R1=${accession}_1.fastq.gz
    R2=${accession}_2.fastq.gz

    api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    report=\$(curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors "\$api" | awk -F'\\t' 'NR==2') || true
    urls=\$(echo "\$report" | cut -f2)
    md5s=\$(echo "\$report" | cut -f3)
    url1=\$(echo "\$urls" | cut -d';' -f1)
    url2=\$(echo "\$urls" | cut -d';' -f2)
    md5_1=\$(echo "\$md5s" | cut -d';' -f1)
    md5_2=\$(echo "\$md5s" | cut -d';' -f2)
    if [ -z "\$url1" ] || [ -z "\$url2" ] || [ "\$url1" = "\$url2" ]; then
        echo "ERROR: ENA lists no paired FASTQ files for ${accession}: \$report" >&2
        exit 1
    fi

    # Report the aggregate rate of k parallel 128 MB range requests against ENA.
    probe() {
        k=\$1
        t0=\$(date +%s.%N)
        for i in \$(seq 0 \$((k - 1))); do
            from=\$((i * 134217728)); to=\$((from + 134217727))
            curl -fsSL -r \$from-\$to -o /dev/null "https://\$url2" &
        done
        wait
        t1=\$(date +%s.%N)
        awk -v k=\$k -v a=\$t0 -v b=\$t1 'BEGIN{printf "[FETCH] ENA rate probe: %d connection(s), %.1f MB/s\\n", k, k * 128 / (b - a)}' >&2
    }

    count_ok() {
        got=\$(gzip -dc "\$1" 2>/dev/null | wc -l | tr -d ' ') || got=0
        [ "\$got" -eq ${n_lines} ]
    }

    stream_subset() {
        curl -fsSL --retry 3 --retry-delay 5 --retry-all-errors "https://\$1" \\
            | gzip -dc | head -n ${n_lines} | gzip -1 -c > "\$2" || true
        count_ok "\$2"
    }

    fetch_full() {
        t0=\$(date +%s)
        curl -fsSL --retry 10 --retry-delay 15 --retry-all-errors -o "\$2" "https://\$1" || return 1
        t1=\$(date +%s)
        bytes=\$(stat -c %s "\$2")
        awk -v f="\$2" -v n=\$bytes -v s=\$((t1 - t0)) 'BEGIN{printf "[FETCH] %s: %.1f GB in %d s, %.1f MB/s\\n", f, n / 1e9, s, n / 1e6 / (s > 0 ? s : 1)}' >&2
        echo "\$3  \$2" | md5sum -c - >&2
    }

    if [ ${n_lines} -gt 0 ]; then
        probe 1 || true
        probe 8 || true
        stream_subset "\$url1" "\$R1" && stream_subset "\$url2" "\$R2" \\
            || { echo "ERROR: could not stream ${max_pairs} read pairs of ${accession} from ENA" >&2; exit 1; }
    else
        fetch_full "\$url1" "\$R1" "\$md5_1" &
        pid1=\$!
        fetch_full "\$url2" "\$R2" "\$md5_2" &
        pid2=\$!
        wait \$pid1 || { echo "ERROR: download of \$R1 failed" >&2; exit 1; }
        wait \$pid2 || { echo "ERROR: download of \$R2 failed" >&2; exit 1; }
    fi
    ls -l "\$R1" "\$R2" >&2
    """

    stub:
    """
    printf '@r#1_2_3/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${accession}_1.fastq.gz
    printf '@r#1_2_3/2\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${accession}_2.fastq.gz
    """
}
