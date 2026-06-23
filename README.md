# RNA seq Pipeline and Downstream Analysis: CAF Reprogramming (GSE280564)

A complete RNA seq workflow spanning a custom Nextflow DSL2 pipeline built on Docker and GCP Batch, through DESeq2 differential expression analysis. The project investigates transcriptional reprogramming of normal human fibroblasts (HUFs) into cancer associated fibroblasts (CAFs) in ovarian cancer, comparing in vitro conditioned media models against genuine tumor derived CAFs.

**Live analysis report:** https://rpubs.com/jriya186/caf-resister-genes-ovarian-cancer

## Dataset

* **GEO Accession:** GSE280564
* **Samples:** 12 total across four conditions, 3 replicates each
  * HUF (normal uterine fibroblasts)
  * Kuramochi conditioned CAF
  * SKOV3 conditioned CAF
  * Primary OC CAF (patient derived)
* **Sequencing:** Illumina NovaSeq 6000, paired end, cDNA
* **Biology:** CAF reprogramming in the ovarian cancer microenvironment
* **Reference paper:** Axemaker et al. (2024), *BBA Molecular Cell Research*

## Project Overview

This project has two parts:

1. **Upstream pipeline** (this repo, Nextflow DSL2) — raw FASTQ to gene count matrix
2. **Downstream analysis** (`analysis/` folder, R/DESeq2) — differential expression, pathway enrichment, and a novel resister gene analysis identifying transcriptional changes that conditioned media models fail to recapitulate from genuine tumor derived CAFs

## Pipeline Steps

1. **SRA_FETCH** — downloads FASTQs directly from SRA inside the pipeline (GCP mode only)
2. **FastQC** — raw read quality control
3. **Fastp** — adapter trimming and quality filtering
4. **STAR** — splice aware alignment to full GRCh38
5. **featureCounts** — gene level read counting across all 12 samples simultaneously, producing a single consistent count matrix
6. **MultiQC** — aggregated QC report

## Tools and Containers

All tools run inside Docker containers from `quay.io/biocontainers`, never Docker Hub, with exact version tags pinned for reproducibility.

| Tool | Version | Container |
|------|---------|-----------|
| SRA Tools | 3.0.3 | `quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0` |
| FastQC | 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| Fastp | 0.23.4 | `quay.io/biocontainers/fastp:0.23.4--h5f740d0_0` |
| STAR | 2.7.11a | `quay.io/biocontainers/star:2.7.11a--h0033a41_0` |
| Samtools | 1.19 | `quay.io/biocontainers/samtools:1.19--h50ea8bc_0` |
| featureCounts | 2.0.6 | `quay.io/biocontainers/subread:2.0.6--he4a0461_0` |
| MultiQC | 1.21 | `quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0` |

## Project Structure

```
rnaseq pipeline/
├── main.nf                          workflow logic, local and GCP modes
├── nextflow.config                  local and GCP Batch profiles
├── samplesheet.csv                  sample to condition mapping for GCP mode
├── modules/
│   ├── sra_fetch.nf
│   ├── fastqc.nf
│   ├── fastp.nf
│   ├── star.nf
│   ├── featurecounts.nf
│   └── multiqc.nf
├── scripts/
│   ├── build_star_index.sh          chr22 index for local testing
│   └── build_star_index_gcp.sh      full GRCh38 index built on a GCP VM
└── analysis/
    └── caf_deseq2_analysis.Rmd      DESeq2, pathway enrichment, resister analysis
```

## Usage

### Local testing (chr22 index, subsampled reads)

```bash
nextflow run main.nf -profile local
```

### Resume after failure

```bash
nextflow run main.nf -profile local -resume
```

### Full run on GCP (all 12 samples, full GRCh38, fetched directly from SRA)

```bash
nextflow run main.nf -profile gcp --mode gcp -resume
```

## Requirements

### Local

* Nextflow 25+
* Docker Desktop (with `--platform linux/amd64` for Apple Silicon)
* Java 17+

### GCP

* Google Cloud project with Batch, Compute Engine, and Storage APIs enabled
* GCS bucket for reference files, work directory, and results
* Full GRCh38 STAR index built once on a GCP VM and stored in GCS
* Default compute service account with `roles/batch.agentReporter` and `roles/storage.objectAdmin`

## Reference

* Genome: GRCh38 (Ensembl release 113)
* Annotation: `Homo_sapiens.GRCh38.113.gtf`
* Local testing: chr22 only, 1M subsampled reads per sample via seqtk
* Production: full genome, full read depth, 90%+ mapping rate across all 12 samples

## Downstream Analysis Highlights

The full analysis, including code, plots, and interpretation, is published at https://rpubs.com/jriya186/caf-resister-genes-ovarian-cancer. Key findings:

* A core CAF transcriptional signature of 864 genes is shared across all three CAF types relative to HUF, independently replicating the findings of Axemaker et al. (2024)
* Primary (tumor derived) CAFs show nearly twice as many differentially expressed genes as either conditioned CAF type, reflecting more extensive reprogramming by the full tumor microenvironment
* A novel resister gene analysis identifies over 2,500 genes that primary CAFs change but conditioned media models fail to recapitulate, enriched for developmental positional identity programs (HOX genes, regionalization, pattern specification), suggesting these programs may be epigenetically protected from soluble factor driven reprogramming

## Author Notes

Built as a portfolio project to demonstrate an end to end bioinformatics workflow: dataset selection, containerized pipeline development, cloud infrastructure, and independent statistical analysis with novel biological interpretation beyond the original publication.

## Acknowledgment

This project performs an independent analysis of publicly available RNA seq data deposited alongside the following publication. All credit for the experimental design, sample generation, and original biological characterization belongs to the original authors.

Axemaker H, Plesselova S, Calar K, Jorgensen M, Wollman J, de la Puente P. Reprogramming of normal fibroblasts into ovarian cancer associated fibroblasts via non vesicular paracrine signaling induces an activated fibroblast phenotype. *BBA Molecular Cell Research.* 2024;1871:119801. doi:10.1016/j.bbamcr.2024.119801
