configfile: "config/config.yaml"


rule all:
    input:
        "results/20_annotation/3_pharokka",
        "results/20_annotation/3_pharokka_plots",
        "results/20_annotation/4_diamond_vog/prodigal-gv_vs_vog_annotated.tsv",
        "results/20_annotation/5_terl_tree/TerL_tree.treefile",
        "results/20_annotation/5_terl_tree/itol/itol_vog_group.txt",
        "results/20_annotation/5_terl_tree/itol/itol_lca.txt"


# Step 21a: Extract all-contig cluster representatives >= 1000 bp
rule extract_long_contigs:
    input:
        fasta="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta"
    output:
        fasta="results/20_annotation/1_selet_seq/long_contigs.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/20_annotation/1_selet_seq/extract_long_contigs.log",
        err="log/20_annotation/1_selet_seq/extract_long_contigs.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/20_annotation/1_selet_seq
        mkdir -p log/20_annotation/1_selet_seq
        seqkit seq -m 1000 {input.fasta} > {output.fasta} \
            2> {log.err}
        """

# Step 21b: Extract all-contig cluster representatives detected by any primer (strict or relaxed)
rule extract_primer_detected_contigs:
    input:
        fasta="results/14_all_contig_cluster/3_cluster_fasta/all_contigs_cluster_representatives.fasta",
        strict="results/17_primer_investigation/all_contigs/all_contigs_primersearch_strict.txt",
        relaxed="results/17_primer_investigation/all_contigs/all_contigs_primersearch_relaxed.txt"
    output:
        list="results/20_annotation/1_selet_seq/primer_detected_contigs.txt",
        fasta="results/20_annotation/1_selet_seq/primer_detected_contigs.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/20_annotation/1_selet_seq/extract_primer_detected_contigs.log",
        err="log/20_annotation/1_selet_seq/extract_primer_detected_contigs.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/20_annotation/1_selet_seq
        mkdir -p log/20_annotation/1_selet_seq

        grep -h "Sequence:" {input.strict} {input.relaxed} \
            | awk '{{print $2}}' | sort -u > {output.list} \
            2> {log.err}

        seqkit grep -f {output.list} {input.fasta} > {output.fasta} \
            2>> {log.err}
        """

# Step 21c: Merge and deduplicate the length- and primer-selected contig sets
rule merge_annotation_input_contigs:
    input:
        long_contigs="results/20_annotation/1_selet_seq/long_contigs.fasta",
        primer_contigs="results/20_annotation/1_selet_seq/primer_detected_contigs.fasta"
    output:
        fasta="results/20_annotation/1_selet_seq/annotation_input.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/20_annotation/1_selet_seq/merge_annotation_input_contigs.log",
        err="log/20_annotation/1_selet_seq/merge_annotation_input_contigs.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/20_annotation/1_selet_seq
        cat {input.long_contigs} {input.primer_contigs} | seqkit rmdup -n -o {output.fasta} \
            > {log.out} 2> {log.err}
        """

# Step 21e: Pharokka gene annotation
rule pharokka_annotation:
    input:
        representatives="results/20_annotation/1_selet_seq/annotation_input.fasta"
    output:
        directory("results/20_annotation/3_pharokka")
    conda: "envs/pharokka.yaml"
    log:
        out="log/20_annotation/3_pharokka/pharokka_annotation.log",
        err="log/20_annotation/3_pharokka/pharokka_annotation.err"
    threads: config["pharokka"]["threads"]
    resources:
        slurm_partition=config["pharokka"]["partition"],
        runtime=config["pharokka"]["runtime"],
        mem_mb_per_cpu=config["pharokka"]["mem_mb_per_cpu"],
        cpus_per_task=config["pharokka"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        db=config["pharokka_db"],
        gene_predictor=config["pharokka"]["gene_predictor"],
        meta="--meta" if config["pharokka"]["meta"] else "",
        meta_hmm="--meta_hmm" if config["pharokka"]["meta_hmm"] else ""
    shell:
        """
        mkdir -p log/20_annotation/3_pharokka
        pharokka run \
            -i {input.representatives} \
            -o {output} \
            -d {params.db} \
            -t {threads} \
            {params.meta} \
            -g {params.gene_predictor} \
            {params.meta_hmm} \
            > {log.out} 2> {log.err}
        """

# Step 21e2: Remove GenBank records that crash pharokka multiplot (IndexError in
# np.quantile on an empty length array). Empirically this hits short/sparse contigs
# even when they have a non-hypothetical CDS, so contigs are kept only if they are
# both >= min_contig_length and have >=1 non-hypothetical CDS. A crash on any one
# contig kills the whole multiplot run and drops plots for every later contig too.
rule filter_pharokka_gbk_zero_cds:
    input:
        genbank="results/20_annotation/3_pharokka/pharokka.gbk"
    output:
        genbank="results/20_annotation/3_pharokka_gbk_filtered/pharokka_nonzero_cds.gbk"
    conda: "envs/python.yaml"
    log:
        out="log/20_annotation/3_pharokka_plots/filter_pharokka_gbk_zero_cds.log",
        err="log/20_annotation/3_pharokka_plots/filter_pharokka_gbk_zero_cds.err"
    resources:
        slurm_partition=config["pharokka"]["multiplot"]["partition"],
        runtime=config["pharokka"]["multiplot"]["runtime"],
        mem_mb_per_cpu=config["pharokka"]["multiplot"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    params:
        min_length=config["pharokka"]["multiplot"]["min_contig_length"]
    shell:
        """
        mkdir -p results/20_annotation/3_pharokka_gbk_filtered
        mkdir -p log/20_annotation/3_pharokka_plots
        python workflow/scripts/3_annotation/filter_gbk_zero_cds.py \
            --genbank {input.genbank} \
            --output {output.genbank} \
            --min-length {params.min_length} \
            > {log.out} 2> {log.err}
        """

# Step 21e3: Plot per-contig genome maps from the Pharokka annotation (pharokka multiplot,
# formerly pharokka_multiplotter.py)
# Step 21e3a: Split the filtered GenBank into N chunks for parallel plotting (see below)
rule split_pharokka_gbk_for_multiplot:
    input:
        genbank="results/20_annotation/3_pharokka_gbk_filtered/pharokka_nonzero_cds.gbk"
    output:
        chunk_dir=directory("results/20_annotation/3_pharokka_plots_chunks")
    conda: "envs/python.yaml"
    log:
        out="log/20_annotation/3_pharokka_plots/split_pharokka_gbk_for_multiplot.log",
        err="log/20_annotation/3_pharokka_plots/split_pharokka_gbk_for_multiplot.err"
    resources:
        slurm_partition=config["pharokka"]["multiplot"]["partition"],
        runtime=config["pharokka"]["multiplot"]["runtime"],
        mem_mb_per_cpu=config["pharokka"]["multiplot"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    params:
        chunks=config["pharokka"]["multiplot"]["threads"]
    shell:
        """
        mkdir -p log/20_annotation/3_pharokka_plots
        python workflow/scripts/3_annotation/split_gbk_chunks.py \
            --genbank {input.genbank} \
            --chunks {params.chunks} \
            --outdir {output.chunk_dir} \
            > {log.out} 2> {log.err}
        """

# Step 21e3b: Plot per-contig genome maps from the Pharokka annotation (pharokka multiplot,
# formerly pharokka_multiplotter.py). pharokka multiplot has no --threads option and plots
# one GenBank file's records sequentially in a single process, so parallelism is done here
# by backgrounding one plain process per pre-split chunk within this rule's own CPU
# allocation (not srun sub-steps: on this cluster, nested `srun --exclusive` job-step
# creation hangs indefinitely with "Requested nodes are busy" when other jobs share the
# node, since it tries to claim whole-node resources rather than just this job's slice).
# Each chunk directory holds one record per .gbk file (see split_gbk_chunks.py), and
# multiplot is invoked once per record with `|| true`, so a np.quantile crash on any
# single contig only drops that one PNG instead of killing every later contig in its chunk.
rule pharokka_multiplot:
    input:
        chunk_dir="results/20_annotation/3_pharokka_plots_chunks"
    output:
        directory("results/20_annotation/3_pharokka_plots")
    conda: "envs/pharokka.yaml"
    log:
        out="log/20_annotation/3_pharokka_plots/pharokka_multiplot.log",
        err="log/20_annotation/3_pharokka_plots/pharokka_multiplot.err"
    threads: config["pharokka"]["multiplot"]["threads"]
    resources:
        slurm_partition=config["pharokka"]["multiplot"]["partition"],
        runtime=config["pharokka"]["multiplot"]["runtime"],
        mem_mb_per_cpu=config["pharokka"]["multiplot"]["mem_mb_per_cpu"],
        cpus_per_task=config["pharokka"]["multiplot"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        label_hypotheticals="--label_hypotheticals" if config["pharokka"]["multiplot"]["label_hypotheticals"] else "",
        remove_other_features_labels="--remove_other_features_labels" if config["pharokka"]["multiplot"]["remove_other_features_labels"] else "",
        title_size=config["pharokka"]["multiplot"]["title_size"],
        label_size=config["pharokka"]["multiplot"]["label_size"],
        interval=config["pharokka"]["multiplot"]["interval"],
        truncate=config["pharokka"]["multiplot"]["truncate"],
        dpi=config["pharokka"]["multiplot"]["dpi"],
        annotations=config["pharokka"]["multiplot"]["annotations"]
    shell:
        """
        mkdir -p log/20_annotation/3_pharokka_plots
        mkdir -p {output}

        for chunk_dir in {input.chunk_dir}/chunk_*/; do
            (
                for record in "$chunk_dir"record_*.gbk; do
                    [ -e "$record" ] || continue
                    record_name=$(basename "$record" .gbk)
                    rm -rf "${{chunk_dir}}plots_${{record_name}}"
                    pharokka multiplot \
                        -g "$record" \
                        -o "${{chunk_dir}}plots_${{record_name}}" \
                        -f \
                        {params.label_hypotheticals} \
                        {params.remove_other_features_labels} \
                        --title_size {params.title_size} \
                        --label_size {params.label_size} \
                        --interval {params.interval} \
                        --truncate {params.truncate} \
                        --dpi {params.dpi} \
                        --annotations {params.annotations} \
                        >> {log.out} 2>> {log.err} || echo "SKIPPED (crashed): $record" >> {log.err}
                done
            ) &
        done
        wait

        cp {input.chunk_dir}/chunk_*/plots_*/*.png {output}/ 2>> {log.err} || true
        """

# Step 21f: Build DIAMOND database from VOG protein sequences
rule build_vog_diamond_db:
    input:
        fasta=config["vog"]["db_dir"] + "/vogdb.proteins.all.fa"
    output:
        db=config["vog"]["db_dir"] + "/vogdb.proteins.all.dmnd"
    conda: "envs/diamond.yaml"
    log:
        out="log/20_annotation/4_diamond_vog/build_vog_diamond_db.log",
        err="log/20_annotation/4_diamond_vog/build_vog_diamond_db.err"
    threads: config["diamond_vog"]["threads"]
    resources:
        slurm_partition=config["diamond_vog"]["partition"],
        runtime=config["diamond_vog"]["runtime"],
        mem_mb_per_cpu=config["diamond_vog"]["mem_mb_per_cpu"],
        cpus_per_task=config["diamond_vog"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/20_annotation/4_diamond_vog
        diamond makedb --in {input.fasta} -d {output.db} -p {threads} \
            > {log.out} 2> {log.err}
        """

# Step 21g: DIAMOND BLASTP of Pharokka-predicted proteins (prodigal-gv) against VOG
rule diamond_vog_annotation:
    input:
        query="results/20_annotation/3_pharokka/prodigal-gv.faa",
        db=config["vog"]["db_dir"] + "/vogdb.proteins.all.dmnd"
    output:
        tsv="results/20_annotation/4_diamond_vog/prodigal-gv_vs_vog.tsv"
    conda: "envs/diamond.yaml"
    log:
        out="log/20_annotation/4_diamond_vog/diamond_vog_annotation.log",
        err="log/20_annotation/4_diamond_vog/diamond_vog_annotation.err"
    threads: config["diamond_vog"]["threads"]
    resources:
        slurm_partition=config["diamond_vog"]["partition"],
        runtime=config["diamond_vog"]["runtime"],
        mem_mb_per_cpu=config["diamond_vog"]["mem_mb_per_cpu"],
        cpus_per_task=config["diamond_vog"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        evalue=config["diamond_vog"]["evalue"],
        min_identity=config["diamond_vog"]["min_identity"],
        min_query_cover=config["diamond_vog"]["min_query_cover"]
    shell:
        """
        mkdir -p results/20_annotation/4_diamond_vog
        mkdir -p log/20_annotation/4_diamond_vog
        diamond blastp \
            --query {input.query} \
            --db {input.db} \
            --evalue {params.evalue} \
            --id {params.min_identity} \
            --query-cover {params.min_query_cover} \
            --threads {threads} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
            --out {output.tsv} \
            > {log.out} 2> {log.err}
        """

# Step 21h: Annotate best-hit DIAMOND/VOG matches with VOG group, function, and LCA taxonomy
rule annotate_vog_hits:
    input:
        diamond="results/20_annotation/4_diamond_vog/prodigal-gv_vs_vog.tsv",
        members=config["vog"]["db_dir"] + "/vog.members.tsv",
        annotations=config["vog"]["db_dir"] + "/vog.annotations.tsv",
        lca=config["vog"]["db_dir"] + "/vog.lca.tsv"
    output:
        tsv="results/20_annotation/4_diamond_vog/prodigal-gv_vs_vog_annotated.tsv"
    conda: "envs/python.yaml"
    log:
        out="log/20_annotation/4_diamond_vog/annotate_vog_hits.log",
        err="log/20_annotation/4_diamond_vog/annotate_vog_hits.err"
    resources:
        slurm_partition=config["diamond_vog"]["partition"],
        runtime=config["diamond_vog"]["runtime"],
        mem_mb_per_cpu=config["diamond_vog"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/20_annotation/4_diamond_vog
        python workflow/scripts/3_annotation/annotate_vog_hits.py \
            --diamond {input.diamond} \
            --members {input.members} \
            --annotations {input.annotations} \
            --lca {input.lca} \
            --output {output.tsv} \
            > {log.out} 2> {log.err}
        """

# Step 21i: Select VOG groups annotated as terminase large subunit
rule select_terl_vog_groups:
    input:
        annotations=config["vog"]["db_dir"] + "/vog.annotations.tsv"
    output:
        list="results/20_annotation/5_terl_tree/terl_vog_groups.txt"
    conda: "envs/python.yaml"
    log:
        out="log/20_annotation/5_terl_tree/select_terl_vog_groups.log",
        err="log/20_annotation/5_terl_tree/select_terl_vog_groups.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    params:
        keyword=config["tree_build"]["terl_keyword"]
    shell:
        """
        mkdir -p results/20_annotation/5_terl_tree
        mkdir -p log/20_annotation/5_terl_tree
        python workflow/scripts/3_annotation/select_terl_vogs.py \
            --annotations {input.annotations} \
            --keyword "{params.keyword}" \
            --output {output.list} \
            > {log.out} 2> {log.err}
        """

# Step 21j: Extract per-group VOG reference sequences for the selected TerL VOGs
rule extract_terl_vog_sequences:
    input:
        list="results/20_annotation/5_terl_tree/terl_vog_groups.txt"
    output:
        fasta="results/20_annotation/5_terl_tree/terl_vog_references.faa"
    log:
        out="log/20_annotation/5_terl_tree/extract_terl_vog_sequences.log",
        err="log/20_annotation/5_terl_tree/extract_terl_vog_sequences.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    params:
        faa_dir=config["vog"]["db_dir"] + "/faa"
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        > {output.fasta}
        while read -r vog; do
            cat {params.faa_dir}/${{vog}}.faa >> {output.fasta}
        done < {input.list} \
            > {log.out} 2> {log.err}
        """

# Step 21k: Combine query TerL sequences (from Pharokka) with VOG reference sequences
rule combine_terl_sequences:
    input:
        query="results/20_annotation/3_pharokka/terL.faa",
        references="results/20_annotation/5_terl_tree/terl_vog_references.faa"
    output:
        fasta="results/20_annotation/5_terl_tree/all_TerL.faa"
    log:
        out="log/20_annotation/5_terl_tree/combine_terl_sequences.log",
        err="log/20_annotation/5_terl_tree/combine_terl_sequences.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        cat {input.query} {input.references} > {output.fasta} \
            2> {log.err}
        """

# Step 21k2: Strip trailing stop codons and remove fragment sequences before alignment
# (short fragments and internal-stop artifacts otherwise force MAFFT into huge gap
# blocks that trimAl then strips entirely, producing an empty alignment)
rule filter_terl_sequences:
    input:
        fasta="results/20_annotation/5_terl_tree/all_TerL.faa"
    output:
        fasta="results/20_annotation/5_terl_tree/all_TerL_filtered.faa"
    conda: "envs/seqkit.yaml"
    log:
        out="log/20_annotation/5_terl_tree/filter_terl_sequences.log",
        err="log/20_annotation/5_terl_tree/filter_terl_sequences.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    params:
        min_length=config["tree_build"]["min_seq_length"]
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        seqkit seq -m {params.min_length} {input.fasta} \
            | sed '/^[^>]/ s/\\*$//' > {output.fasta} \
            2> {log.err}
        """

# Step 21l: Align combined TerL sequences with MAFFT (L-INS-i)
rule align_terl_mafft:
    input:
        fasta="results/20_annotation/5_terl_tree/all_TerL_filtered.faa"
    output:
        fasta="results/20_annotation/5_terl_tree/aligned_TerL.faa"
    conda: "envs/tree_build.yaml"
    log:
        out="log/20_annotation/5_terl_tree/align_terl_mafft.log",
        err="log/20_annotation/5_terl_tree/align_terl_mafft.err"
    threads: config["tree_build"]["threads"]
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        cpus_per_task=config["tree_build"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        maxiterate=config["tree_build"]["mafft_maxiterate"]
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        mafft \
            --localpair \
            --maxiterate {params.maxiterate} \
            --thread {threads} \
            {input.fasta} > {output.fasta} \
            2> {log.err}
        """

# Step 21m: Trim the TerL alignment with trimAl
rule trim_terl_alignment:
    input:
        fasta="results/20_annotation/5_terl_tree/aligned_TerL.faa"
    output:
        fasta="results/20_annotation/5_terl_tree/trimmed_TerL.faa"
    conda: "envs/tree_build.yaml"
    log:
        out="log/20_annotation/5_terl_tree/trim_terl_alignment.log",
        err="log/20_annotation/5_terl_tree/trim_terl_alignment.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        trimal \
            -in {input.fasta} \
            -out {output.fasta} \
            -automated1 \
            > {log.out} 2> {log.err}
        """

# Step 21n: Build the TerL phylogenetic tree with IQ-TREE2
rule build_terl_tree:
    input:
        fasta="results/20_annotation/5_terl_tree/trimmed_TerL.faa"
    output:
        treefile="results/20_annotation/5_terl_tree/TerL_tree.treefile"
    conda: "envs/tree_build.yaml"
    log:
        out="log/20_annotation/5_terl_tree/build_terl_tree.log",
        err="log/20_annotation/5_terl_tree/build_terl_tree.err"
    threads: config["tree_build"]["threads"]
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        cpus_per_task=config["tree_build"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        model=config["tree_build"]["iqtree_model"],
        bootstrap=config["tree_build"]["iqtree_bootstrap"],
        prefix="results/20_annotation/5_terl_tree/TerL_tree"
    shell:
        """
        mkdir -p log/20_annotation/5_terl_tree
        iqtree2 \
            -s {input.fasta} \
            -m {params.model} \
            -B {params.bootstrap} \
            -T AUTO \
            --prefix {params.prefix} \
            > {log.out} 2> {log.err}
        """

# Step 21o: Build iTOL annotation files for the TerL tree's VOG reference leaves
# (VOG group membership and LCA taxonomy, one TREE_COLORS file each)
rule build_terl_itol_annotations:
    input:
        references="results/20_annotation/5_terl_tree/terl_vog_references.faa",
        alignment_input="results/20_annotation/5_terl_tree/all_TerL_filtered.faa",
        members=config["vog"]["db_dir"] + "/vog.members.tsv",
        lca=config["vog"]["db_dir"] + "/vog.lca.tsv"
    output:
        vog_colors="results/20_annotation/5_terl_tree/itol/itol_vog_group.txt",
        lca_colors="results/20_annotation/5_terl_tree/itol/itol_lca.txt"
    conda: "envs/python.yaml"
    log:
        out="log/20_annotation/5_terl_tree/build_terl_itol_annotations.log",
        err="log/20_annotation/5_terl_tree/build_terl_itol_annotations.err"
    resources:
        slurm_partition=config["tree_build"]["partition"],
        runtime=config["tree_build"]["runtime"],
        mem_mb_per_cpu=config["tree_build"]["mem_mb_per_cpu"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/20_annotation/5_terl_tree/itol
        mkdir -p log/20_annotation/5_terl_tree
        python workflow/scripts/3_annotation/build_terl_itol_annotations.py \
            --references {input.references} \
            --alignment-input {input.alignment_input} \
            --members {input.members} \
            --lca {input.lca} \
            --vog-output {output.vog_colors} \
            --lca-output {output.lca_colors} \
            > {log.out} 2> {log.err}
        """
