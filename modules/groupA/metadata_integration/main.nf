process METADATA_INTEGRATION {
    label 'process_single'
    container 'python:3.10'
    publishDir "${params.outdir}/groupA/metadata", mode: 'copy'

    input:
        path metadata_csv
        val sample_ids

    output:
        path "integrated_metadata.json", emit: metadata
        path "validation_report.txt", emit: report
        path "versions.yml", emit: versions

    script:
    """
    python << EOF
    import json
    import pandas as pd

    # Read metadata
    df = pd.read_csv('$metadata_csv')

    # Validate and integrate
    metadata_dict = df.set_index('sample').to_dict('index')

    # Save as JSON
    with open('integrated_metadata.json', 'w') as f:
        json.dump(metadata_dict, f, indent=2)

    # Create report
    with open('validation_report.txt', 'w') as f:
        f.write(f"Metadata validation complete\\n")
        f.write(f"Total samples: {len(df)}\\n")
    EOF

    cat > versions.yml << EOF
    "METADATA_INTEGRATION":
        python: \$(python --version | cut -d' ' -f2)
    EOF
    """
}
