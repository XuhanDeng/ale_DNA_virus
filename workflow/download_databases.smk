configfile: "config/config.yaml"


rule all:
    input:
        config["checkv_db"] + "done.txt",
        config["genomad_db"] + "_done.txt",
        config["virsorter2"]["database"] + "/done.txt",
        config["pharokka_db"] + "/done.txt",
        config["vog"]["db_dir"] + "/done.txt"


# Database 1: CheckV
rule download_checkv_db:
    output:
        marker=config["checkv_db"] + "done.txt"
    conda: "envs/checkv.yaml"
    log:
        out="log/download_databases/checkv.log",
        err="log/download_databases/checkv.err"
    threads: config["checkv"]["threads"]
    shell:
        """
        mkdir -p database
        mkdir -p log/download_databases

        checkv download_database database \
            > {log.out} 2> {log.err}

        touch {output.marker}
        """


# Database 2: geNomad
rule download_genomad_db:
    output:
        marker=config["genomad_db"] + "_done.txt"
    conda: "envs/genomad.yaml"
    log:
        out="log/download_databases/genomad.log",
        err="log/download_databases/genomad.err"
    threads: config["genomad"]["threads"]
    params:
        db_parent=lambda wc: config["genomad_db"].rsplit("/", 1)[0]
    shell:
        """
        mkdir -p {params.db_parent}
        mkdir -p log/download_databases

        genomad download-database {params.db_parent} \
            > {log.out} 2> {log.err}

        touch {output.marker}
        """


# Database 3: VirSorter2
rule download_virsorter2_db:
    output:
        marker=config["virsorter2"]["database"] + "/done.txt"
    conda: "envs/vs2.yaml"
    log:
        out="log/download_databases/virsorter2.log",
        err="log/download_databases/virsorter2.err"
    threads: config["virsorter2"]["threads"]
    params:
        db_dir=config["virsorter2"]["database"]
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p log/download_databases

        virsorter setup -d {params.db_dir} -j {threads} \
            > {log.out} 2> {log.err}

        touch {output.marker}
        """


# Database 4: Pharokka
# Uses pharokka's own installer so the database version always matches the
# installed pharokka version (avoids the version-mismatch that manual
# Zenodo tarball downloads are prone to as pharokka releases new versions).
rule download_pharokka_db:
    output:
        marker=config["pharokka_db"] + "/done.txt"
    conda: "envs/pharokka.yaml"
    log:
        out="log/download_databases/pharokka.log",
        err="log/download_databases/pharokka.err"
    params:
        db_dir=config["pharokka_db"]
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p log/download_databases

        pharokka install -o {params.db_dir} \
            > {log.out} 2> {log.err}

        touch {output.marker}
        """


# Database 5: VOG (Virus Orthologous Groups)
rule download_vog_db:
    output:
        marker=config["vog"]["db_dir"] + "/done.txt"
    log:
        out="log/download_databases/vog.log",
        err="log/download_databases/vog.err"
    params:
        db_dir=config["vog"]["db_dir"],
        version=config["vog"]["version"]
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p log/download_databases

        wget -c -P {params.db_dir} https://fileshare.csb.univie.ac.at/vog/{params.version}/vog.lca.tsv.gz \
            > {log.out} 2> {log.err}
        wget -c -P {params.db_dir} https://fileshare.csb.univie.ac.at/vog/{params.version}/vog.faa.tar.gz \
            >> {log.out} 2>> {log.err}
        wget -c -P {params.db_dir} https://fileshare.csb.univie.ac.at/vog/{params.version}/vog.annotations.tsv.gz \
            >> {log.out} 2>> {log.err}
        wget -c -P {params.db_dir} https://fileshare.csb.univie.ac.at/vog/{params.version}/vog.members.tsv.gz \
            >> {log.out} 2>> {log.err}
        wget -c -P {params.db_dir} https://fileshare.csb.univie.ac.at/vog/{params.version}/vogdb.proteins.all.fa.gz \
            >> {log.out} 2>> {log.err}

        gunzip -f {params.db_dir}/vog.lca.tsv.gz \
            >> {log.out} 2>> {log.err}
        gunzip -f {params.db_dir}/vog.annotations.tsv.gz \
            >> {log.out} 2>> {log.err}
        gunzip -f {params.db_dir}/vog.members.tsv.gz \
            >> {log.out} 2>> {log.err}
        gunzip -f {params.db_dir}/vogdb.proteins.all.fa.gz \
            >> {log.out} 2>> {log.err}
        tar -xzf {params.db_dir}/vog.faa.tar.gz -C {params.db_dir} \
            >> {log.out} 2>> {log.err}

        rm -f {params.db_dir}/vog.faa.tar.gz

        touch {output.marker}
        """
