# 1. Load the required software onto the compute node
module purge
module load nextflow
module load singularityCE/4.3

# 2. Apptainer / Singularity image cache
export NXF_APPTAINER_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/
export NXF_SINGULARITY_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/

# 3. Run
echo "=== Launching nf-core/abgtemplate ==="
nextflow run main.nf \
    -profile groupB,slurm,singularity \
    --input data/sample_sheet.csv \
    --metadata data/patients_sample_information.tsv \
    --outdir results_integ \
    -ansi-log \
    -resume

# Demo smoke uses the bundled 5-row data/sample_sheet.csv (raw .fastq under data/seqs/,
# (schema: sample,fastq_1[,fastq_2]) built from the fasterq-dump output, e.g.: