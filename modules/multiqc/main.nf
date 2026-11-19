// Collect the per-sample statistics into one MultiQC report. The pipeline's own JSON and
// text statistics are first turned into MultiQC custom-content tables by a module binary;
// samtools markdup, flagstat and stats files are parsed by MultiQC itself.
process MULTIQC {
    conda 'bioconda::multiqc=1.25'

    input:
    path qc_files
    path config

    output:
    path 'multiqc_report.html', emit: report
    // MultiQC names the data directory after --filename.
    path 'multiqc_report_data', emit: data

    script:
    """
    stlfr_mqc_tables.py .
    multiqc --force --config ${config} --filename multiqc_report.html .
    """

    stub:
    """
    mkdir -p multiqc_report_data
    echo '<html><body>stub</body></html>' > multiqc_report.html
    """
}
