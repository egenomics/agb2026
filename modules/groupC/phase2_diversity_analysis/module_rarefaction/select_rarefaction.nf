process RAREFACTION_THRESHOLD {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.11 numpy pandas matplotlib"

    input:
    tuple val(meta), path(asv_table)

    output:
    tuple val(meta), path("${meta.id}/rarefaction_threshold.txt"), emit: threshold
    tuple val(meta), path("${meta.id}/sample_qc.tsv"), emit: qc
    tuple val(meta), path("${meta.id}/report.txt"), emit: report
    tuple val(meta), path("${meta.id}/*.pdf"), emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def method = task.ext.method ?: 'knee'
    def percentile = task.ext.percentile ?: 10
    def cov_target = task.ext.coverage_target ?: 0.99
    def dropout_max = task.ext.dropout_max ?: 0.10
    def step = task.ext.step ?: 500
    def iterations = task.ext.iterations ?: 100

    """
    mkdir -p ${meta.id}

    pip install --quiet --no-cache-dir numpy pandas matplotlib

    python ${projectDir}/select_rarefaction_depth.py \\
        --input ${asv_table} \\
        --output ${meta.id} \\
        --method ${method} \\
        --percentile ${percentile} \\
        --coverage-target ${cov_target} \\
        --dropout-max ${dropout_max} \\
        --step ${step} \\
        --iterations ${iterations} \\
        ${args}
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo "10000" > ${meta.id}/rarefaction_threshold.txt
    echo -e "sample\tlib_size\tpasses_threshold\tcoverage" > ${meta.id}/sample_qc.tsv
    echo "Stub report" > ${meta.id}/report.txt

    """
}
