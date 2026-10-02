# Freshwater Viral Metagenome Pipeline

A Snakemake pipeline for viral identification, quality filtering, clustering, abundance quantification, and taxonomic annotation from paired-end metagenomic reads. Designed for SLURM HPC clusters with Conda environment management.

## Pipeline Overview

```text
Raw reads (FASTQ)
    │
    ▼ Step 1  fastp — quality control & adapter trimming
    │
    ▼ Step 4  SPAdes (--meta) — metagenome assembly
    │
    ▼ Step 5  seqkit — rename contigs, filter by min length (default 5000 bp)
    │
    ├──▶ GeNomad — viral identification
    └──▶ VirSorter2 — viral identification
    │
    ▼ Step 10  Merge GeNomad + VirSorter2 results & pre-filter
    │
    ▼ Step 11  seqkit — extract candidate viral contigs
    │
    ▼ Step 12  CheckV — genome quality assessment
    │
    ▼ Step 13  Python — second filter based on CheckV results
    │
    ▼ Step 14  seqkit — extract final high-confidence viral contigs
    │
    ├──▶ Cluster per sample (BLAST + ANI)
    └──▶ Cluster all samples combined (BLAST + ANI)  ◀── active path
              │
              ▼ Extract representative sequences (cluster centroids)
              │
              ▼ Build Bowtie2 index
              │
              ▼ Bowtie2 alignment (all samples → shared index)
              │
              ▼ CoverM — abundance: TPM, mean coverage, read count
              │
              ▼ GeNomad annotate — taxonomy for cluster representatives
              │
              ▼ iPHoP — host prediction for cluster representatives
```

> **Note:** Per-sample abundance/annotation paths are implemented but commented out in the current Snakefile. Only the all-samples combined path is active.

## Requirements

- Snakemake ≥ 7.0
- Conda / Mamba
- SLURM cluster

### Tools (managed via Conda environments in `workflow/envs/`)

| Tool | Purpose | Conda env |
| ---- | ------- | --------- |
| fastp | Read QC | `workflow/envs/fastp.yaml` |
| SPAdes | Assembly | `workflow/envs/spades.yaml` |
| seqkit | Sequence filtering/extraction | `workflow/envs/seqkit.yaml` |
| GeNomad | Viral identification & taxonomy | `workflow/envs/genomad.yaml` |
| VirSorter2 | Viral identification | `workflow/envs/vs2.yaml` |
| CheckV | Viral genome quality | `workflow/envs/checkv.yaml` |
| BLAST + Python | ANI-based clustering | `workflow/envs/checkv.yaml` |
| Bowtie2 | Read alignment | `workflow/envs/bowtie2.yaml` |
| Samtools | BAM processing | `workflow/envs/bowtie2_samtools.yaml` |
| CoverM | Abundance calculation | `workflow/envs/coverm.yaml` |
| iPHoP | Host prediction | `workflow/envs/iphop.yaml` |

## Input

Paired-end gzipped FASTQ files placed in `input/`:

```text
input/
├── {sample}.R1.raw.fastq.gz
└── {sample}.R2.raw.fastq.gz
```

## Configuration

All parameters are set in `config/config.yaml`. Key settings to update before running:

```yaml
# List your sample names
samples:
  - sample1
  - sample2

# SLURM
slurm_account: "your-account"

# Database paths
checkv_db: "database/checkv-db-v1.5/"
genomad_db: "path/to/genomad_db"
virsorter2:
  database: "path/to/virsorter2_db"

# Script paths
merge_script: "workflow/scripts/merge_viral_identification.py"
checkv_filter_script: "workflow/scripts/Second_filter.py"
checkv_scripts: "workflow/scripts/checkv"
```

### Key filtering parameters

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `seqkit.min_length` | 5000 | Minimum contig length (bp) after assembly |
| `virsorter2.min_score` | 0.5 | VirSorter2 minimum score threshold |
| `merge_filter.vs2_min_score` | 0.5 | Pre-filter score cutoff when merging tools |
| `virsorter2.groups` | dsDNAphage, NCLDV, lavidaviridae | Target viral groups |

## iPHoP Database Setup

Before running host prediction, download and set up the iPHoP database using `workflow/iphop_establish.smk`. This is a one-time setup step.

```bash
snakemake --snakefile workflow/iphop_establish.smk --profile config/slurm/
```

What it does:

1. Creates `database/iphop/` and `database/iphop_test/` directories
2. Downloads `iPHoP.latest_rw` from NERSC portal (the current published version's actual
   directory name, e.g. `Jun_2025_pub_rw`, may differ from this download alias and will
   change whenever NERSC publishes a new "latest" release)
3. Verifies MD5 checksum integrity
4. Extracts the database archive

The downloaded database path should match `iphop.database` in `config/config.yaml`. Check the
end of `log/iphop_establish/download.err` (it reports "The database has been put in ...") to
confirm the actual extracted directory name before setting this:

```yaml
iphop:
  database: "database/iphop/Jun_2025_pub_rw"
```

## Running the Pipeline

```bash
# Dry run (check pipeline graph without executing)
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --profile config/slurm/

# Submit to SLURM (see RUN_COMMANDS.md for the full command)
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile config/slurm/
```

The profile enables: SLURM executor, up to 50 concurrent jobs, Conda environments, and automatic re-run of incomplete jobs.

## Output Structure

```text
results/
├── 1_fastp/           # QC reads
├── 3_spades_result/   # Assembly scaffolds
├── 4_rename_assembly/ # Filtered & renamed contigs
├── 5_genomad/         # GeNomad viral predictions
├── 6_virsorter2/      # VirSorter2 viral predictions
├── 8_python_merge_filter/  # Merged results & candidate viral contigs
├── 9_checkv/          # CheckV quality summaries & filtered lists
├── 10_cluster/
│   ├── 1_filtered_seq/     # Final per-sample viral contigs
│   ├── 2_combined/         # All samples combined
│   ├── 3_clustered/        # BLAST + ANI clustering results
│   ├── 4_cluster_fasta/    # Cluster representative sequences
│   └── 5_bowtie_index/     # Bowtie2 indices
├── 11_bowtie2/
│   └── all_samples/{sample}/
│       ├── {sample}_TPM.tsv      # TPM abundance
│       ├── {sample}_mean.tsv     # Mean coverage
│       └── {sample}_count.tsv   # Read counts
├── 12_taxonomy/genomad/all_samples/  # Taxonomic annotation
└── 13_host_prediction/iphop/         # iPHoP host prediction results
```

Logs are written to `log/` mirroring the results structure.

## Detailed Pipeline Steps

### Step 1 — Quality Control (`fastp_quality_control`)

Raw paired-end reads are quality-trimmed using fastp. Poly-G and poly-X tails are trimmed, and reads failing quality or length thresholds are discarded. Both paired and unpaired outputs are retained. QC reports are generated in HTML and JSON format.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Minimum base quality (Phred) | 20 | `fastp.qualified_quality_phred` |
| Minimum read length | 20 bp | `fastp.length_required` |
| Threads | 8 | `fastp.threads` |

### Step 4 — Metagenome Assembly (`spades_assembly`)

Quality-filtered reads are assembled using SPAdes in metagenome mode (`--meta`). Multiple k-mer sizes are used to maximize recovery of diverse viral genomes. Temporary files and intermediate graph files are removed after assembly to save disk space.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| k-mer sizes | 21,33,55,77,101,121 | `spades.k_values` |
| Max memory | 256 GB | `spades.memory` |
| Threads | 16 | `spades.threads` |

### Step 5 — Rename and Filter Assemblies (`rename_filter_assemblies`)

Assembled scaffolds are filtered by minimum length and renamed with a consistent `{sample}_{number}` scheme using seqkit. This ensures unique contig identifiers across samples for downstream merging and clustering.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Minimum contig length | 5000 bp | `seqkit.min_length` |
| ID number width | 10 digits | `seqkit.nr_width` |

### Steps 6–7 — Viral Identification (GeNomad & VirSorter2)

Two complementary viral identification tools are run in parallel on the filtered contigs:

- **GeNomad** (`genomad_identification`): uses a marker-based approach combined with a neural network to classify sequences as viral, plasmid, or chromosomal. Run in end-to-end mode with `--cleanup`.
- **VirSorter2** (`virsorter2_identification`): uses hallmark gene detection and machine learning to identify viral sequences.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| VirSorter2 minimum score | 0.5 | `virsorter2.min_score` |
| VirSorter2 minimum length | 5000 bp | `virsorter2.min_length` |
| VirSorter2 target groups | dsDNAphage, NCLDV, lavidaviridae | `virsorter2.groups` |

### Step 10 — Merge and Pre-filter Viral Results (`merge_viral_results`)

Outputs from GeNomad and VirSorter2 are merged using a custom Python script (`workflow/scripts/merge_viral_identification.py`). Sequences are retained if detected by either tool. The merged result table and a list of candidate viral contig IDs are written for the next step.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| VirSorter2 minimum score | 0.5 | `merge_filter.vs2_min_score` |
| VirSorter2 allow partial sequences | false (full only) | `merge_filter.vs2_allow_partial` |
| GeNomad exclude topology | NONE (all topologies kept) | `merge_filter.genomad_exclude_topology` |

`vs2_allow_partial: false` means only VirSorter2 hits flagged as `full` are kept. Set to `true` to also include partial detections.
`genomad_exclude_topology` accepts a topology string to exclude (e.g. `"Provirus"`). Set to `"NONE"` to keep all topologies.

### Step 11 — Extract Candidate Viral Sequences (`extract_viral_sequences`)

seqkit greps the candidate viral contig IDs from Step 10 out of the renamed assembly FASTA. This produces a reduced FASTA to pass to CheckV for quality assessment.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Minimum contig length | 5000 bp | `seqkit.min_length` |
| Threads | 4 | `seqkit.threads` |

### Step 12 — CheckV Quality Assessment (`checkv_filtered_analysis`)

CheckV evaluates the completeness and quality of each candidate viral contig, identifying complete, high-quality, medium-quality, and low-quality genomes. It also trims host contamination from provirus sequences. Results are written to `quality_summary.tsv`.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Database | checkv-db-v1.5 | `checkv_db` |
| Threads | 4 | `checkv.threads` |

### Step 13 — Second Filter (`python_checkv_filter`)

A second round of filtering is applied using `workflow/scripts/Second_filter.py`, integrating CheckV quality results with the original GeNomad and VirSorter2 scores. Each tool's results are filtered independently using the criteria below, then the passing sequences are union-merged. This removes low-confidence predictions and produces per-tool filtered tables plus a final extraction list.

| Filter | Criterion | Description |
| ------ | --------- | ----------- |
| VirSorter2 hallmark | `hallmark > 0` | Requires at least one hallmark viral gene |
| VirSorter2 completeness | `full/partial == "full"` | Retains only full (non-partial) sequences |
| GeNomad topology | `topology != "Provirus"` | Excludes provirus predictions |
| GeNomad hallmark | `n_hallmarks > 0` | Requires at least one hallmark gene |
| CheckV viral genes | `viral_genes > 0` | Requires at least one viral gene called by CheckV |

These thresholds are hard-coded in `workflow/scripts/Second_filter.py` and are not controlled by `config/config.yaml`.

**Retention logic is OR (union):** a sequence is kept if it passes the criteria for **any one** of the three tools — it does NOT need to satisfy all three. For example, a contig with `viral_genes > 0` in CheckV alone is sufficient for retention, even if it was not detected by VirSorter2 or GeNomad.

**Outputs** (written to `results/9_checkv/{sample}/`):

| File | Content |
| ---- | ------- |
| `{sample}_checkv_extract.txt` | Final list of contig IDs to extract (union of all three filters) |
| `{sample}_vs2_filtered.tsv` | VirSorter2 sequences passing the second filter |
| `{sample}_genomad_filtered.tsv` | GeNomad sequences passing the second filter |
| `{sample}_checkv_filtered.tsv` | CheckV sequences passing the second filter |
| `{sample}_source_table.tsv` | Per-contig table indicating which tool(s) passed each sequence |

### Step 14 — Extract Final Viral Sequences (`extract_final_viral_sequences`)

The final high-confidence viral contigs are extracted from the renamed assembly FASTA using the list from Step 13. These per-sample FASTA files are the input to all downstream clustering and abundance steps.

### Steps 15a/b — ANI-based Clustering (`cluster_per_sample` / `cluster_all_samples`)

Viral contigs are dereplicated into species-level clusters (viral OTUs) using the CheckV ANI-based method:

1. All-vs-all BLASTn search
2. ANI calculation with `anicalc.py`
3. Greedy clustering with `aniclust.py`

Both per-sample and cross-sample (all samples combined) clustering are run. Only the all-samples path feeds into downstream steps.

| Parameter | Value |
| --------- | ----- |
| Minimum ANI | 95% |
| Minimum target coverage | 85% |
| Minimum query coverage | 0% |
| Threads | 16 |

### Step 16b — Extract Cluster Representatives (`extract_representatives_all_samples`)

The centroid sequence from each ANI cluster is extracted using seqkit. These representative sequences define the non-redundant viral population catalogue used for all downstream steps.

### Step 17b — Build Bowtie2 Index (`build_cluster_index_all_samples`)

A Bowtie2 alignment index is built from the cluster representative sequences. All samples are mapped against this shared index to enable cross-sample abundance comparison.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Threads | 4 | `cluster_index_build.threads` |

### Step 18b — Read Alignment (`bowtie2_alignment_all_samples`)

Quality-filtered reads from each sample are aligned to the shared Bowtie2 index in sensitive mode. Only properly paired reads (`-f 2`) are retained, converted to sorted BAM, and indexed. SAM files are deleted after processing to save disk space.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Alignment mode | `--sensitive` | — |
| Retained reads | Properly paired only (`-f 2`) | — |
| Threads | 8 | `bowtie2_alignment.threads` |

### Step 19b — Abundance Calculation (`calculate_abundance_all_samples`)

CoverM calculates three abundance metrics per contig per sample from the sorted BAM files.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Minimum read identity | 95% | — |
| Minimum aligned fraction | 90% | — |
| Outputs | TPM, mean coverage, read count | — |
| Threads | 4 | `coverm.threads` |

### Step 20b — Taxonomic Annotation (`genomad_annotate_all_samples`)

GeNomad is run in `annotate` mode on the cluster representative sequences to assign viral taxonomy at family and genus level based on the geNomad database.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Database | geNomad DB | `genomad_db` |
| Threads | 16 | `genomad.threads` |

### Step 21 — Host Prediction (`iphop_host_prediction`)

Run separately via `workflow/3_annotation_host_prediction.smk` (see `RUN_COMMANDS.md`), since it only depends on the viral cluster representatives and has no downstream consumers in the main pipeline.

iPHoP predicts the bacterial or archaeal host for each viral cluster representative by integrating CRISPR spacer matching, alignment to reference genomes, and sequence composition. Results are written to `results/13_host_prediction/iphop/`.

| Parameter | Value | Config key |
| --------- | ----- | ---------- |
| Database | Jun_2025_pub_rw (extracted from `iPHoP.latest_rw` download) | `iphop.database` |
| Threads | 32 | `iphop.threads` |

## Clustering Method

Viral contigs are dereplicated using the CheckV ANI-based approach:

1. All-vs-all BLASTn
2. ANI calculation (`anicalc.py`)
3. Greedy clustering (`aniclust.py`) with:
   - `--min_ani 95`
   - `--min_tcov 85`
   - `--min_qcov 0`

Cluster centroids (representatives) are used for all downstream abundance and taxonomy steps.

## Abundance Calculation

CoverM is run on BAM files with strict alignment filters:

- Minimum read identity: 95%
- Minimum aligned fraction: 90%
- Outputs: TPM, mean coverage, and raw read count per contig

## Citations

If you use this pipeline, please cite:

- **fastp**: Chen et al. (2018) *Bioinformatics*
- **SPAdes**: Bankevich et al. (2012) *J. Comput. Biol.*
- **GeNomad**: Camargo et al. (2023) *Nature Biotechnology*
- **VirSorter2**: Guo et al. (2021) *Microbiome*
- **CheckV**: Nayfach et al. (2021) *Nature Biotechnology*
- **CoverM**: github.com/wwood/CoverM
- **seqkit**: Shen et al. (2016) *PLOS ONE*
- **iPHoP**: Roux et al. (2023) *PLOS Biology*
