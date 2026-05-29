# 1. Load the required software onto the compute node
module purge
module load nextflow
module load singularityCE/4.3

# 2. Apptainer / Singularity image cache
export NXF_APPTAINER_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/
export NXF_SINGULARITY_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/

# 3. Run
echo "=== Llançant Nextflow Test ==="
nextflow run test_groupB.nf \
    -profile groupB,slurm,singularity \
    --input assets/samplesheet.csv \
    --outdir results \
    -ansi-log \
    -resume
