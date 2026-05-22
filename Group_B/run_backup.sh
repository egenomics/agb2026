#!/bin/bash

# 1. Environment Hygiene: Load the required software onto the compute node
module purge
module load nextflow
module load singularityCE/4.3

# 2. Definim la variable de la memòria cau d'Apptainer/Singularity
export NXF_APPTAINER_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/
export NXF_SINGULARITY_CACHEDIR=/data/upfagb/jvillanueva/apptainer_cache/

# 3. Neteja de la memòria cau dels processos que van fallar amb errors o fitxers buits
echo "=== Netejant directoris de treball conflictius ==="
nextflow clean -w 27/2d5a42d54a75a50a5f9 -f
nextflow clean -w fd/fe789412d487af797a7896b945da96 -f

# 4. Execute your isolated Sandbox Test
echo "=== Llançant Nextflow Test ==="
nextflow run test_groupB.nf \
    -profile cluster,singularity \
    --input assets/samplesheet_test.csv \
    --outdir results_test \
    -ansi log \
