configfile: "config/config.yaml"


rule all:
    input:
        "results/99_statistic/merged_abundance/all_samples_TPM.csv",
        "results/99_statistic/merged_abundance/all_samples_TPM.xlsx",
        "results/99_statistic/merged_abundance/all_samples_mean.csv",
        "results/99_statistic/merged_abundance/all_samples_mean.xlsx",
        "results/99_statistic/merged_abundance/all_samples_count.csv",
        "results/99_statistic/merged_abundance/all_samples_count.xlsx"


# Step 1: Compute contig length for all-contig cluster representatives
rule contig_length:
    input:
        representatives="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta"
    output:
        length="results/99_statistic/contig_length.tsv"
    conda: "envs/seqkit.yaml"
    log:
        out="log/99_statistic/contig_length.log",
        err="log/99_statistic/contig_length.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/99_statistic
        mkdir -p log/99_statistic

        echo -e "Contig\\tlength" > {output.length}
        seqkit fx2tab -n -l {input.representatives} >> {output.length} \
            2> {log.err}
        """

# Step 2: Parse primersearch strict/relaxed results into a per-contig primer tag table
rule primer_tags:
    input:
        strict="results/17_primer_investigation/all_contigs/all_contigs_primersearch_strict.txt",
        relaxed="results/17_primer_investigation/all_contigs/all_contigs_primersearch_relaxed.txt"
    output:
        tags="results/99_statistic/primer_tags.tsv"
    conda: "envs/python.yaml"
    log:
        out="log/99_statistic/primer_tags.log",
        err="log/99_statistic/primer_tags.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/99_statistic
        mkdir -p log/99_statistic

        python workflow/scripts/99_statistic/parse_primersearch.py \
            --strict {input.strict} \
            --relaxed {input.relaxed} \
            --output {output.tags} \
            > {log.out} 2> {log.err}
        """

# Step 2b: Parse BLASTN single-primer results into a per-contig single-primer tag table
rule single_primer_tags:
    input:
        blast="results/17_primer_investigation/all_contigs/all_contigs_blast_single_primer.txt"
    output:
        tags="results/99_statistic/single_primer_tags.tsv"
    conda: "envs/python.yaml"
    log:
        out="log/99_statistic/single_primer_tags.log",
        err="log/99_statistic/single_primer_tags.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/99_statistic
        mkdir -p log/99_statistic

        python workflow/scripts/99_statistic/parse_blast_single_primer.py \
            --blast {input.blast} \
            --output {output.tags} \
            > {log.out} 2> {log.err}
        """

# Step 3: Merge per-sample TPM/mean/count abundance tables, annotated with
# length, virus membership, and primer-detection tags; write CSV + XLSX
rule merge_abundance_tables:
    input:
        tpm=expand("results/15_all_contig_bowtie2/{sample}/{sample}_TPM.tsv", sample=config["samples"]),
        mean=expand("results/15_all_contig_bowtie2/{sample}/{sample}_mean.tsv", sample=config["samples"]),
        read_count=expand("results/15_all_contig_bowtie2/{sample}/{sample}_count.tsv", sample=config["samples"]),
        length="results/99_statistic/contig_length.tsv",
        virus_fasta="results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta",
        primer_tags="results/99_statistic/primer_tags.tsv",
        single_primer_tags="results/99_statistic/single_primer_tags.tsv"
    output:
        tpm_csv="results/99_statistic/merged_abundance/all_samples_TPM.csv",
        tpm_xlsx="results/99_statistic/merged_abundance/all_samples_TPM.xlsx",
        mean_csv="results/99_statistic/merged_abundance/all_samples_mean.csv",
        mean_xlsx="results/99_statistic/merged_abundance/all_samples_mean.xlsx",
        count_csv="results/99_statistic/merged_abundance/all_samples_count.csv",
        count_xlsx="results/99_statistic/merged_abundance/all_samples_count.xlsx"
    conda: "envs/python.yaml"
    log:
        out="log/99_statistic/merge_abundance_tables.log",
        err="log/99_statistic/merge_abundance_tables.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    params:
        samples=config["samples"]
    shell:
        """
        mkdir -p results/99_statistic/merged_abundance
        mkdir -p log/99_statistic

        python workflow/scripts/99_statistic/merge_abundance.py \
            --inputs {input.tpm} \
            --samples {params.samples} \
            --length {input.length} \
            --virus-fasta {input.virus_fasta} \
            --primer-tags {input.primer_tags} \
            --single-primer-tags {input.single_primer_tags} \
            --output-csv {output.tpm_csv} \
            --output-xlsx {output.tpm_xlsx} \
            > {log.out} 2> {log.err}

        python workflow/scripts/99_statistic/merge_abundance.py \
            --inputs {input.mean} \
            --samples {params.samples} \
            --length {input.length} \
            --virus-fasta {input.virus_fasta} \
            --primer-tags {input.primer_tags} \
            --single-primer-tags {input.single_primer_tags} \
            --output-csv {output.mean_csv} \
            --output-xlsx {output.mean_xlsx} \
            >> {log.out} 2>> {log.err}

        python workflow/scripts/99_statistic/merge_abundance.py \
            --inputs {input.read_count} \
            --samples {params.samples} \
            --length {input.length} \
            --virus-fasta {input.virus_fasta} \
            --primer-tags {input.primer_tags} \
            --single-primer-tags {input.single_primer_tags} \
            --output-csv {output.count_csv} \
            --output-xlsx {output.count_xlsx} \
            >> {log.out} 2>> {log.err}
        """
