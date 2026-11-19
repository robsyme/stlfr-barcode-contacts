# stlfr-barcode-contacts

A Nextflow pipeline that turns stLFR linked reads into barcode co-occurrence contact maps, one per sample and chromosome.

## Background

stLFR (single-tube long fragment read) libraries label all the short reads that come from one long DNA molecule with the same barcode. If two positions on a chromosome are close together, many molecules span both and the two positions share many barcodes. The number of shared barcodes falls off with distance and reaches a flat background beyond the molecule length, which is usually 100 to 300 kb.

The pipeline aligns the reads to a reference, cuts each chromosome into 5 kb bins, and counts for every pair of bins the barcodes seen in both. Plotted as a matrix, a correct assembly shows a band along the diagonal. Strong signal away from the diagonal means that two regions the reference places far apart sit on the same molecules in the sample, which is what an inversion, translocation, duplication, collapsed repeat or misassembly looks like.

## Steps

1. `FETCH_READS` downloads the paired FASTQ files of an SRA run from ENA. With `--max_pairs N` it streams only the first N read pairs; otherwise it downloads the full files and checks their md5 sums. With `--read_source sra`, `FETCH_SRA` runs instead: it downloads the run's `.sra` object from the location NCBI reports, checks its md5 sum and converts it to FASTQ with one `fastq-dump` per CPU. Use it where ENA is slow to reach.
2. `TAG_BARCODES` reads the 42 bp barcode block at the end of read 2, matches the three 10 bp barcodes against the whitelist, writes them into the read name (`name#A_B_C`), trims read 2 to the 100 bp insert and drops pairs without a valid barcode. It can keep a random fraction of barcodes (`--barcode_fraction`), which thins the data without breaking up molecules. Output is split into chunks of `--chunk_pairs` read pairs.
3. `PREPARE_REFERENCE` and `BWAMEM2_INDEX` decompress and index the reference.
4. `ALIGN` runs `bwa-mem2 mem`, `samtools fixmate -m` and `samtools sort` on each chunk.
5. `MERGE_MARKDUP` merges the chunks of a sample and removes duplicates with `samtools markdup -r`.
6. `BARCODE_MATRIX` runs once per sample and chromosome. It takes alignments with MAPQ of at least 20, records the barcode and bin (`int(POS / bin_size)`) of each, and counts the distinct barcodes shared by every pair of bins.
7. `MAKE_BCM` converts each table into a `.bcm` bundle for the [barcode contact map viewer](https://github.com/robsyme/presentations/tree/main/2026-11-barcode-contact-map), which also holds the hotspot and junction calls.
8. `BAM_STATS` runs `samtools flagstat` and `samtools stats` on each duplicate-removed BAM, and `MULTIQC` gathers these with the tagging, duplicate and table statistics into one report.

## Inputs

A sample sheet in CSV format with a header:

```csv
sample,accession
PBA_Drummond,SRR32381415
PBA_HatTrick,SRR32381414
```

Sample names are used in file names, so keep them free of spaces. The accession is an SRA run accession with paired FASTQ files on ENA.

The reference is a FASTA file, plain, gzipped or bgzipped.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--input` | required | Sample sheet |
| `--reference` | required | Reference FASTA |
| `--chromosomes` | `^Cec[0-9]+$` | Regular expression; tables are built for reference sequences whose whole name matches |
| `--read_source` | `ena` | `ena` or `sra`: where to download the reads from |
| `--max_pairs` | `0` | Read pairs taken from the start of each run; 0 means the full run |
| `--barcode_fraction` | `0.25` | Fraction of barcodes kept |
| `--seed` | `1` | Seed for the barcode subsample |
| `--chunk_pairs` | `10000000` | Read pairs per alignment chunk |
| `--bin_size` | `5000` | Bin size in bp |
| `--barcode_whitelist` | bundled | stLFR barcode whitelist |
| `--publish_bam` | `false` | Also publish the duplicate-removed BAMs |

## Running

The pipeline needs Nextflow 25.10 or later. Each process declares its software with a `conda` directive (one uses a prebuilt container), and the helper scripts are module binaries, so the usual way to run it is with Wave building the containers:

```bash
nextflow run robsyme/stlfr-barcode-contacts \
    -profile wave \
    --input samples.csv \
    --reference reference.fa.gz \
    -output-dir results
```

Nothing in the repository names a bucket, queue or compute environment. Give the work directory, output directory and executor at launch. On Seqera Platform, launch it on a compute environment that has Wave enabled and set `outputDir` in the launch configuration.

For a quick test, take a few million read pairs and keep every barcode:

```bash
nextflow run robsyme/stlfr-barcode-contacts -profile wave \
    --input samples.csv --reference reference.fa.gz \
    --max_pairs 4000000 --barcode_fraction 1.0 --chunk_pairs 1500000 \
    -output-dir results-test
```

`nextflow run . -stub -profile local_small` checks the workflow wiring without running any tools.

## Outputs

Outputs are declared in the workflow `output` block.

```
bcm/<sample>/<chrom>.bcm                              viewer bundle
stat/<sample>/<sample>.<chrom>_5kb_bin.stat.gz        co-occurrence table
stat/<sample>/<sample>.<chrom>_5kb_bin.summary.json   record, barcode and count totals
multiqc/multiqc_report.html                           MultiQC report (data in multiqc/multiqc_report_data/)
bam_stats/<sample>.flagstat, <sample>.stats           samtools statistics of the duplicate-removed BAM
tag_stats/<sample>.tag_stats.json                     barcode tagging statistics
markdup_stats/<sample>.markdup_stats.txt              samtools markdup statistics
bam/<sample>/<sample>.markdup.bam                     only with --publish_bam
```

The table has three tab-separated integer columns without a header: start of bin A, start of bin B, and the number of barcodes found in both. It is symmetric: each pair appears in both orders and the diagonal appears once. The diagonal holds the number of barcodes in that bin.

## Viewing

The viewer is a single static HTML page. It loads bundles from a file picker or from URLs given as `index.html?src=<url>`, with `src` repeated once per bundle. A bundle stores its chromosome name and the viewer keeps one bundle per chromosome name, so load the bundles of one sample at a time. The bundles of different samples have the same file names and are kept apart by their `bcm/<sample>/` directory.

## Data used during development

The pipeline was developed on the stLFR runs of NCBI BioProject PRJNA1225167 (Garg et al. 2025, "An Australian chickpea pan-genome provides insights into genome organization", Plant Biotechnology Journal 23:3967-3983, https://doi.org/10.1111/pbi.70192) aligned to a *Cicer echinospermum* reference whose chromosomes are named `Cec1` to `Cec8`. Read 1 is 100 bp and read 2 is 142 bp: 100 bp of insert followed by the barcode block.

## Tests

`tests/` holds unit tests for the barcode tagging script. Run them with `uv run --with pytest pytest tests/`.
