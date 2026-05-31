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
    --dada2_train_set /data/upfagb/u269668/agb2026/biodb/dada2/silva_nr99_v138.2_toGenus_trainset.fa \
    --dada2_species_set /data/upfagb/u269668/agb2026/biodb/dada2/silva_v138.2_assignSpecies.fa \
    -ansi-log \
    -resume
