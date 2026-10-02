configfile: "config/config.yaml"


rule all:
    input:
        "results/17_primer_investigation/all_contigs/all_contigs_primersearch_strict.txt",
        "results/17_primer_investigation/all_contigs/all_contigs_primersearch_relaxed.txt",
        "results/17_primer_investigation/viral_cluster/viral_cluster_primersearch_strict.txt",
        "results/17_primer_investigation/viral_cluster/viral_cluster_primersearch_relaxed.txt",
        "results/17_primer_investigation/all_contigs/all_contigs_blast_single_primer.txt"


# Step 0: Convert primer FASTA to EMBOSS primersearch -infile format
rule convert_primers_to_primersearch_format:
    input:
        fasta="input/primer.fasta"
    output:
        infile="results/17_primer_investigation/primers.primersearch"
    conda: "envs/python.yaml"
    log:
        out="log/28_primer_investigation/convert_primers_to_primersearch_format.log",
        err="log/28_primer_investigation/convert_primers_to_primersearch_format.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/17_primer_investigation
        mkdir -p log/28_primer_investigation
        python workflow/scripts/fasta_to_primersearch.py \
            --fasta {input.fasta} \
            --output {output.infile} \
            > {log.out} 2> {log.err}
        """

# Step 1a: primersearch against all-contig cluster representatives (strict, 0% mismatch)
rule primersearch_all_contigs_strict:
    input:
        fasta="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta",
        infile="results/17_primer_investigation/primers.primersearch"
    output:
        result="results/17_primer_investigation/all_contigs/all_contigs_primersearch_strict.txt"
    conda: "envs/emboss.yaml"
    log:
        out="log/28_primer_investigation/primersearch_all_contigs_strict.log",
        err="log/28_primer_investigation/primersearch_all_contigs_strict.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/17_primer_investigation/all_contigs
        mkdir -p log/28_primer_investigation
        primersearch -seqall {input.fasta} \
            -infile {input.infile} \
            -mismatchpercent 0 \
            -outfile {output.result} \
            > {log.out} 2> {log.err}
        """

# Step 1b: primersearch against all-contig cluster representatives (relaxed, 10% mismatch)
rule primersearch_all_contigs_relaxed:
    input:
        fasta="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta",
        infile="results/17_primer_investigation/primers.primersearch"
    output:
        result="results/17_primer_investigation/all_contigs/all_contigs_primersearch_relaxed.txt"
    conda: "envs/emboss.yaml"
    log:
        out="log/28_primer_investigation/primersearch_all_contigs_relaxed.log",
        err="log/28_primer_investigation/primersearch_all_contigs_relaxed.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/17_primer_investigation/all_contigs
        mkdir -p log/28_primer_investigation
        primersearch -seqall {input.fasta} \
            -infile {input.infile} \
            -mismatchpercent 10 \
            -outfile {output.result} \
            > {log.out} 2> {log.err}
        """

# Step 1c: build a BLAST nucleotide database from the all-contig cluster representatives
rule build_all_contigs_blast_db:
    input:
        fasta="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta"
    output:
        db="results/17_primer_investigation/all_contigs_blastdb/all_contigs.ndb"
    conda: "envs/blast.yaml"
    log:
        out="log/28_primer_investigation/build_all_contigs_blast_db.log",
        err="log/28_primer_investigation/build_all_contigs_blast_db.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    params:
        db_prefix="results/17_primer_investigation/all_contigs_blastdb/all_contigs"
    shell:
        """
        mkdir -p results/17_primer_investigation/all_contigs_blastdb
        mkdir -p log/28_primer_investigation
        makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_prefix} \
            > {log.out} 2> {log.err}
        """

# Step 1d: single-primer search (BLASTN) against all-contig cluster representatives.
# Unlike primersearch, this does not require a paired F/R primer bounding an amplicon -
# each primer (forward or reverse) is searched on its own, in both orientations
# (BLASTN searches both strands by default), reporting every contig a single primer
# alone can hit. task=blastn-short is BLAST's recommended mode for short (<50bp) queries
# like primers. Strict/relaxed tiering is applied downstream from pident/mismatch columns.
rule blast_single_primer_all_contigs:
    input:
        primers="input/primer.fasta",
        db="results/17_primer_investigation/all_contigs_blastdb/all_contigs.ndb"
    output:
        result="results/17_primer_investigation/all_contigs/all_contigs_blast_single_primer.txt"
    conda: "envs/blast.yaml"
    log:
        out="log/28_primer_investigation/blast_single_primer_all_contigs.log",
        err="log/28_primer_investigation/blast_single_primer_all_contigs.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    params:
        db_prefix="results/17_primer_investigation/all_contigs_blastdb/all_contigs"
    shell:
        """
        mkdir -p results/17_primer_investigation/all_contigs
        mkdir -p log/28_primer_investigation
        blastn -task blastn-short \
            -query {input.primers} \
            -db {params.db_prefix} \
            -word_size 7 \
            -evalue 1000 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
            -out {output.result} \
            > {log.out} 2> {log.err}
        """

# Step 2a: primersearch against viral cluster representatives (strict, 0% mismatch)
rule primersearch_viral_cluster_strict:
    input:
        fasta="results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta",
        infile="results/17_primer_investigation/primers.primersearch"
    output:
        result="results/17_primer_investigation/viral_cluster/viral_cluster_primersearch_strict.txt"
    conda: "envs/emboss.yaml"
    log:
        out="log/28_primer_investigation/primersearch_viral_cluster_strict.log",
        err="log/28_primer_investigation/primersearch_viral_cluster_strict.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/17_primer_investigation/viral_cluster
        mkdir -p log/28_primer_investigation
        primersearch -seqall {input.fasta} \
            -infile {input.infile} \
            -mismatchpercent 0 \
            -outfile {output.result} \
            > {log.out} 2> {log.err}
        """

# Step 2b: primersearch against viral cluster representatives (relaxed, 10% mismatch)
rule primersearch_viral_cluster_relaxed:
    input:
        fasta="results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta",
        infile="results/17_primer_investigation/primers.primersearch"
    output:
        result="results/17_primer_investigation/viral_cluster/viral_cluster_primersearch_relaxed.txt"
    conda: "envs/emboss.yaml"
    log:
        out="log/28_primer_investigation/primersearch_viral_cluster_relaxed.log",
        err="log/28_primer_investigation/primersearch_viral_cluster_relaxed.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/17_primer_investigation/viral_cluster
        mkdir -p log/28_primer_investigation
        primersearch -seqall {input.fasta} \
            -infile {input.infile} \
            -mismatchpercent 10 \
            -outfile {output.result} \
            > {log.out} 2> {log.err}
        """
