# Run all commands from the project working directory (e.g. /scratch/xunwang/aleida)

# ===================== Setup (once per shell session) =====================
XDG_CACHE_HOME="$(pwd)/.cache"
export XDG_CACHE_HOME
mkdir -p "$XDG_CACHE_HOME"


# ===================== workflow/download_databases.smk =====================
# Run on LOGIN NODE with internet access (no slurm executor)
# NOTE: this also installs iPHoP from source (workflow/scripts/install_iphop_from_source.sh),
# which creates its own "iphop" conda env (via git clone + mamba env create, not
# Snakemake's conda: auto-provisioning). workflow/iphop_establish.smk and
# workflow/4_iphop.smk both `conda activate iphop` directly in their shell blocks and
# require this env to already exist, so run download_databases.smk to completion first.

# dry-run
snakemake -s workflow/download_databases.smk --use-conda --cores 1 --rerun-triggers input -n -p

# build conda envs only
snakemake -s workflow/download_databases.smk --sdm conda --conda-create-envs-only

# run (login node)
mkdir -p log/download_databases
nohup snakemake -s workflow/download_databases.smk \
    --use-conda \
    --rerun-triggers input \
    --cores 4 \
    --latency-wait 10 \
    > log/download_databases/snakemake.log 2>&1 &
echo "PID: $!"

# unlock
snakemake -s workflow/download_databases.smk --unlock


# ===================== workflow/iphop_establish.smk =====================
# Run on LOGIN NODE with internet access (no slurm executor)

# dry-run
snakemake -s workflow/iphop_establish.smk --use-conda --cores 1 --rerun-triggers input -n -p

# build conda envs only
snakemake -s workflow/iphop_establish.smk --sdm conda --conda-create-envs-only

# run (login node)
mkdir -p log/iphop_establish
nohup snakemake -s workflow/iphop_establish.smk \
    --use-conda \
    --rerun-triggers input \
    --cores 1 \
    --latency-wait 10 \
    > log/iphop_establish/snakemake.log 2>&1 &
echo "PID: $!"

# unlock
snakemake -s workflow/iphop_establish.smk --unlock


# ===================== workflow/Snakefile =====================

# dry-run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 1 -n -p --rerun-triggers input

# build conda envs only
snakemake -s workflow/Snakefile --configfile config/config.yaml --sdm conda --conda-create-envs-only

# unlock
snakemake -s workflow/Snakefile --configfile config/config.yaml --unlock

# formal run
mkdir -p log/Snakefile
nohup snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --executor slurm \
    --jobs 48 \
    --use-conda \
    --retries 2 \
    --printshellcmds \
    --slurm-no-account \
    --keep-going \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/Snakefile/snakemake.log 2>&1 &
echo "PID: $!"


# ===================== workflow/2_primer_investigation.smk =====================

# dry-run
snakemake -s workflow/2_primer_investigation.smk --configfile config/config.yaml --use-conda --cores 1 -n -p --rerun-triggers input

# build conda envs only
snakemake -s workflow/2_primer_investigation.smk --configfile config/config.yaml --sdm conda --conda-create-envs-only

# unlock
snakemake -s workflow/2_primer_investigation.smk --configfile config/config.yaml --unlock

# formal run
mkdir -p log/28_primer_investigation
nohup snakemake -s workflow/2_primer_investigation.smk --configfile config/config.yaml \
    --executor slurm \
    --jobs 4 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 10 \
    > log/28_primer_investigation/snakemake.log 2>&1 &
echo "PID: $!"


# ===================== workflow/3_annotation.smk =====================

# dry-run
snakemake -s workflow/3_annotation.smk --configfile config/config.yaml --use-conda --cores 1 -n -p --rerun-triggers input 

# build conda envs only
snakemake -s workflow/3_annotation.smk --configfile config/config.yaml --sdm conda --conda-create-envs-only

# unlock
snakemake -s workflow/3_annotation.smk --configfile config/config.yaml --unlock

# formal run
mkdir -p log/20_annotation
nohup snakemake -s workflow/3_annotation.smk --configfile config/config.yaml \
    --executor slurm \
    --jobs 4 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 10 \
    > log/20_annotation/3_annotation_snakemake.log 2>&1 &
echo "PID: $!"


# ===================== workflow/4_iphop.smk =====================

# dry-run
snakemake -s workflow/4_iphop.smk --configfile config/config.yaml --use-conda --cores 1 -n -p --rerun-triggers input

# build conda envs only
snakemake -s workflow/4_iphop.smk --configfile config/config.yaml --sdm conda --conda-create-envs-only

# unlock
snakemake -s workflow/4_iphop.smk --configfile config/config.yaml --unlock

# formal run
mkdir -p log/20_annotation
nohup snakemake -s workflow/4_iphop.smk --configfile config/config.yaml \
    --executor slurm \
    --jobs 4 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 10 \
    > log/20_annotation/snakemake.log 2>&1 &
echo "PID: $!"


# ===================== workflow/99_statistic.smk =====================

# dry-run
snakemake -s workflow/99_statistic.smk --configfile config/config.yaml --use-conda --cores 2 --rerun-triggers input -n -p 

# build conda envs only
snakemake -s workflow/99_statistic.smk --configfile config/config.yaml --sdm conda --conda-create-envs-only

# unlock
snakemake -s workflow/99_statistic.smk --configfile config/config.yaml --unlock

# formal run
mkdir -p log/99_statistic
nohup snakemake -s workflow/99_statistic.smk --configfile config/config.yaml \
    --executor slurm \
    --jobs 4 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 10 \
    > log/99_statistic/snakemake.log 2>&1 &
echo "PID: $!"
