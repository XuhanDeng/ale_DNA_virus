configfile: "config/config.yaml"


rule all:
    input:
        "results/20_annotation/2_iphop"


# Step 21d: iPHoP host prediction
rule iphop_host_prediction:
    input:
        representatives="results/20_annotation/1_selet_seq/annotation_input.fasta"
    output:
        directory("results/20_annotation/2_iphop")
    conda: "envs/iphop.yaml"
    log:
        out="log/20_annotation/2_iphop/iphop_host_prediction.log",
        err="log/20_annotation/2_iphop/iphop_host_prediction.err"
    threads: config["iphop"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["iphop"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        db=config["iphop"]["database"]
    shell:
        """
        mkdir -p results/20_annotation
        mkdir -p log/20_annotation/2_iphop
        iphop predict --fa_file {input.representatives} --out_dir {output} --db_dir {params.db} -t {threads} \
            --max_thread_wish {threads} \
            > {log.out} 2> {log.err}
        """
