configfile: "config/config.yaml"

IPHOP_DB_VERSION = config["iphop"]["db_version"]
IPHOP_DB_DOWNLOAD_ALIAS = config["iphop"]["db_download_alias"]


rule all:
    input:
        expand("database/iphop/" + IPHOP_DB_VERSION + "/{file}", file=["metadata.tsv", "genome_info.tsv", "host_info.tsv", "host_prediction.tsv", "host_prediction_summary.tsv"])


#iphop

rule host_identification_iphop_database_download:
    output:
        expand("database/iphop/" + IPHOP_DB_VERSION + "/{file}", file=["metadata.tsv", "genome_info.tsv", "host_info.tsv", "host_prediction.tsv", "host_prediction_summary.tsv"])
    conda: "envs/iphop.yaml"
    log:
        out="log/iphop_establish/download.log",
        err="log/iphop_establish/download.err"
    shell:
        """
        mkdir -p database/iphop
        mkdir -p log/iphop_establish

        iphop download --no_prompt --split --db_dir database/iphop -dbv {IPHOP_DB_DOWNLOAD_ALIAS} \
            > {log.out} 2> {log.err}

        # Extraction succeeded (outputs above exist) — remove the downloaded archive
        # chunks and checksum to free disk space; only the extracted database is needed
        rm -f database/iphop/{IPHOP_DB_DOWNLOAD_ALIAS}.tar.gz* \
            >> {log.out} 2>> {log.err}
        """
