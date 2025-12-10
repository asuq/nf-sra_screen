# nf-sra_screen


[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.8-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
<!-- [![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/) -->

## Introduction

**nf-sra_screen** is a Nextflow pipeline for taxon‑focused screening and assembly of public SRA runs and local FASTQ files, followed by taxonomical annotation and binning.

![nf-sra_screen diagram](./images/nf-sra_screen_diagram.png)

Given:
- a list of SRA accessions and/or local FASTQ files,
- a pinned NCBI taxonomy snapshot and NCBI <-> GTDB mapping tables,
- a UniProt DIAMOND database
- (Optional) Sandpiper and SingleM marker‑gene databases,

the pipeline will:
1. Discover and filter appropriate SRR runs from SRA metadata (SRA mode).
2. (Optional) pre‑screen samples using Sandpiper and/or SingleM against a GTDB‑derived phylum list.
3. Assemble only read sets that pass the taxonomic screen (short‑read, ONT, PacBio CLR and HiFi).
4. Annotate contigs with DIAMOND and BlobToolKit.
5. (Optional) extract contigs matching user‑specified taxa into per‑taxon FASTA and ID lists.
6. (Optional) Run multiple binners (MetaBAT2, ComeBin, SemiBin, Rosella) and integrate them with DAS Tool.
7. Collate a per‑sample `summary.tsv` with counts and rich failure/success notes, and post‑annotate it using scheduler info from the Nextflow `trace.tsv`.

You can use the pipeline in four modes:
- **Assembly only**: just give SRA/FASTQ + `--taxdump` + `--uniprot_db`.
- **Assembly + binning** (`--binning`): just give SRA/FASTQ + `--taxdump` + `--uniprot_db`.
- **Assembly + taxon screening** (`--taxa`): additionally provide GTDB mapping, and SingleM/Sandpiper databases.
- **Assembly + taxon screening + binning** (`--taxa` & `--binning`): additionally provide GTDB mapping, and SingleM/Sandpiper databases.

The **top‑level orchestration** is split into four named workflows in `main.nf`:
- `PRE_SCREENING` – SRA metadata -> SRR selection -> optional Sandpiper/SingleM screening.
- `ASSEMBLY` – Assembly, DIAMOND, BlobToolKit, optional taxon extraction.
- `BINNING` – MetaBAT2, ComeBin, SemiBin, Rosella, DAS Tool, and binning note aggregation.
- `SUMMARY` – merges all success and failure notes into the final global `summary.tsv`.

## Installation

### Software requirements

- **Nextflow**: `== 25.04.8`
- **Plugins**:
  - `nf-boost@~0.6.0` (used for intermediate files clean‑up and helper functions such as `groupKey`)
- **Container back‑end**:
  - Docker, or
  - Singularity / Apptainer


### Database requirements

- `--taxdump`        NCBI taxdump dir (`nodes.dmp`, `names.dmp`, `taxidlineage.dmp` or classical taxdump)
- `--uniprot_db`     UniProt DIAMOND database (`.dmnd`) (See the [BlobToolKit documentation](https://blobtoolkit.genomehubs.org/install/) for how to build this)
- `--gtdb_ncbi_map`  (with `--taxa`) Dir with NCBI -> GTDB crosswalk: `ncbi_vs_gtdb_bacteria.xlsx`, `ncbi_vs_gtdb_archaea.xlsx`, `gtdb_r226.dic` from [GTDB download](https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/auxillary_files/)
- `--sandpiper_db`   (with `--taxa`) Sandpiper db with `sandpiper_sra.txt`, `sandpiper1.0.0.condensed.tsv`
- `--singlem_db`     (with `--taxa`) SingleM metapackage (e.g. `S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb`)
All tools used by the pipeline are provided via containers defined in nextflow.config.


## Usage

The pipeline can ingest SRA accessions and/or local FASTQ files in the same run.
Internally these are merged before assembly.

### 1. Input: SRA samplesheet (`--sra`)

Prepare a CSV with a single column `sra`, each row representing an SRA project (or study‑level) accession:

`sra.csv`
```{csv}
sra
PRJNAXXXXXX
SRPXXXXXX
ERRXXXXXX
```

Each row can be a project, study, or run
The pipeline will query metadata and expand each project into multiple SRR runs internally.

The metadata files are written under:
```
<outdir>/metadata/<sra>
```

### 2. Input: FASTQ list (`--fastq`)

Prepare a TSV describing local FASTQ files.

`fastq.tsv`
```{tsv}
sample	read_type	reads
A98	hifi	a98.fastq.gz
B27	short	read_1.fastq.gz,read_2.fastq.gz
C03	nanopore	c03_pass.fastq.gz
D48	pacbio	d48.fastq.gz
```

- `sample`: logical sample identifier.
- `read_type`: types of reads used for assembler selection. Use descriptive labels below
    - `short`: short pair-end reads (Illumina, BGISEQ, DNBSEQ) (metaSPAdes),
    - `nanopore`: Nanopore reads (metaFlye),
    - `pacbio`: PacBio CLR reads (metaFlye),
    - `hifi`: PacBio HiFi reads (myloasm).
- `reads`: comma‑separated list of FASTQ paths (absolute or relative); at least one file per row is required. Two or more files are treated as paired‑end for SingleM/metaSPAdes, one as single‑end.


In FASTQ + screening mode (`--fastq_tsv` + `--taxa`), each sample is treated as
```{txt}
sra = sample
srr = sample
platform = UNKNOWN
model = read_type
strategy = UNKNOWN
assembler = read_type
```

SingleM is run in place of Sandpiper for these samples.

In FASTQ + no‑screening mode (`--fastq_tsv`), reads go straight into assembly and optionally binning.


### 3. Input: taxa list (`--taxa`)

Provide a CSV of target taxa if you want taxon‑specific screening and contig extraction:

`taxa.csv`
```{csv}
rank,taxa
phylum,Bacillota
class,Gammaproteobacteria
order,o__Chloroflexales
genus,g__Escherichia
```

Allowed ranks (case-insensitive)
```{csv}
realm,domain,superkingdom,kingdom,phylum,class,order,family,genus,species
```

> [!IMPORTANT]
> - If you do not supply `--taxa`, the pipeline skips SingleM/Sandpiper and taxon‑specific extraction.
> - If you supply the taxon in GTDB style, the pipeline runs SingleM/Sandpiper but skips taxon‑specific extraction


### Full example command
```{bash}
nextflow run asuq/nf-sra_screen \
  -profile <docker/singularity/local/slurm/...> \
  --sra sra.csv \
  --fastq_tsv fastq.tsv \
  --taxdump /path/to/ncbi_taxdump_dir \
  --uniprot_db /path/to/uniprot.dmnd \
  --taxa taxa.csv \
  --gtdb_ncbi_map /path/to/ncbi_vs_gtdb_xlsx_dir \
  --sandpiper_db /path/to/sandpiper_db_dir \
  --singlem_db /path/to/singlem_metapackage \
  --outdir nf-sra_screen_results
```

### Key parameters
- `-profile`         nextflow profile (see below)
- `--sra`            CSV with column `sra` listing project accessions
- `--fastq_tsv`      TSV with columns (`sample,read_type,reads`) listing sample reads
- `--taxdump`        Directory containing NCBI taxdump files; `jsonify_taxdump.py` will create `taxdump.json`
- `--uniprot_db`     UniProt DIAMOND database (`.dmnd`) (Follow [blobtools tutorial](https://blobtoolkit.genomehubs.org/install/))
- `--taxa`           (Optional) CSV with rank,taxa (NCBI or GTDB names). Use it if you want taxonomy screening
- `--gtdb_ncbi_map`  (Optional) Directory with ncbi_vs_gtdb_bacteria.xlsx and ncbi_vs_gtdb_archaea.xlsx. For taxonomy screening
- `--sandpiper_db`   (Optional) Directory with Sandpiper summary tables. For taxonomy screening
- `--singlem_db`     (Optional) SingleM metapackage (e.g. S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb) For taxonomy screening
- `--outdir`         Output directory (default: ./output)
- `--max_retries`    Maximum number of retries per process (default: 3)
- `--help`           Print the pipeline help message and exit.

### Profiles
- `local`
  - Executor: `local`
  - `docker.enabled = true`
  - Small queue size and moderate resources (max_cpus=8, max_memory=16.GB).
- `slurm`
  - Executor: `slurm`
  - `singularity.enabled = true`
  - Large queue size (`queueSize=2000`) and increased resource caps.
- `oist`
  - Includes `conf/oist.config` for OIST Deigo HPC settings.
- `debug`
  - docker.enabled = true
  - `executor.queueSize = 1`
  - Extended trace.fields for debugging.
- `test`
  - For small regression tests.


## Output structure
```{txt}
<output>/
  metadata/
    <sra>/
      <sra>.filtered.csv
      <sra>.skipped.csv
      <sra>.FAIL.note              # if metadata step failed

  <sra>/<srr>/
    # Screening
    sandpiper_report.txt
    sandpiper_output.tsv
    sandpiper_decision.txt
    singlem_taxonomic_profile.tsv
    singlem_taxonomic_profile_krona*
    singlem_output.tsv

    # Assembly
    assembly.fasta
    assembly.gfa
    spades.log / flye.log / myloasm.log
    fastp.html                     # short-read only

    # BlobToolKit
    blobtools.csv
    blobtools*.svg

    # Taxon extraction (if --taxa)
    summary.csv
    *.ids.csv
    *.fasta

    # Binning (if --binning)
    binning/
      metabat/
      comebin/
      semibin/
      rosella/
      dastool/
      metabat.note                 # if failed
      comebin.note                 # if failed
      semibin.note                 # if failed
      rosella.note                 # if failed
      dastool.note                 # if failed
      binning_note.txt             # aggregated notes

  summary.tsv                      # global summary across all samples
  execution-reports/
    timeline.html
    report.html
    trace.tsv
```

## Credits

`Author / maintainer`: Akito Shima (ASUQ), akito-shima[at]oist.jp

### Core tools (via containers)
- iSeq
- SRA toolkit
- Sandpiper
- SingleM
- DIAMOND
- BlobToolKit
- fastp
- metaSPAdes / SPAdes
- Flye
- myloasm
- bowtie2
- minimap2
- samtools
- MetaBAT2
- ComeBin
- SemiBin
- Rosella
- DAS Tool

<!-- ## Citations -->

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use ASUQ/busco_phylogenomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->
