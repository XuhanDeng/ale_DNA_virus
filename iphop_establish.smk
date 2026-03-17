configfile: "profiles/config.yaml"


rule all:
    input:  
expand("database/iphop_test/iPHoP_db_rw_1.4_for-test/{file}", file=["metadata.tsv", "genome_info.tsv", "host_info.tsv", "host_prediction.tsv", "host_prediction_summary.tsv"])


#iphop

rule host_identification_iphop_database_download:
    mkdir -p  database/{iphop,iphop_test}
    iphop download --db_dir database/iphop -dbv iPHoP_db_Jun25_rw

    wget -c -P database/iphop https://portal.nersc.gov/cfs/m342/iphop/db/iPHoP_db_Jun25_rw.tar.gz

    wget -c -P database/iphop https://portal.nersc.gov/cfs/m342/iphop/db/iPHoP_db_Jun25_rw.tar.gz.md5

    # Verify integrity using md5 checksum and extract only if successful
    cd database/iphop
    if md5sum -c iPHoP_db_Jun25_rw.tar.gz.md5; then
        echo "MD5 checksum passed. Extracting database..."
        tar -zxvf iPHoP_db_Jun25_rw.tar.gz
    else
        echo "ERROR: MD5 checksum failed! Database download is corrupted." >&2
        exit 1
    fi
    cd ../..
