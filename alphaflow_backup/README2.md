# Organization of the folders

````
AGB2026/
├── alphaflow/                <-- Group B: The Nextflow Pipeline
│   ├── main.nf               <-- Entry point
│   ├── nextflow.config       <-- Configurations (Profiles, Params)
│   ├── modules/              <-- Local & nf-core modules
│   ├── workflows/            <-- Your alphaflow.nf logic
│   ├── subworkflows/         <-- Helper subworkflows
│   ├── assets/
│   │   ├── samplesheet.csv   <-- Group A fills this in
│   │   └── schema_input.json <-- Validation rules
│   └── .gitignore            <-- MUST include "biodb/" and "work/"
│
├── group_a_scripts/          <-- Group A: Pre-processing scripts
│   ├── trimming.sh
│   └── quality_control.py
│
├── docs/                     <-- Documentation, images for the report
│   └── design_v1.png
│
└── README.md                 <-- General instructions for the whole class
````

## Considerations for Group A

Crucial Point: Group A should NOT upload the .fastq.gz files to GitHub. GitHub is for code, not for large biological data. If they try to upload gigabytes of DNA sequences, the upload will fail or the repository will become unusable.

What Group A should do instead:

- Place the data locally: They keep the FASTQ files on their own computers (or the cluster) in a folder like raw_data/.

- Update the Samplesheet: They edit the assets/samplesheet.csv inside your alphaflow folder on GitHub.

- The paths: In that CSV, they will put the paths to where the files live on the Cluster (since that’s where you’ll likely run the final analysis).

## Code only for group B

This checks if you have any syntax errors or missing files in your Nextflow logic:

````
nf-core pipelines lint
````


You added new variables to nextflow.config (the database paths), but the "security guard" file (nextflow_schema.json) doesn't recognize them.

The Fix:
You need to tell the schema that these parameters exist. You can do this automatically:

````
# Fix DADA2 Denoising
curl -L https://raw.githubusercontent.com/nf-core/ampliseq/2.9.0/modules/local/dada2_denoising.nf -o modules/local/dada2/dada2_denoising.nf

# Fix DADA2 Taxonomy
curl -L https://raw.githubusercontent.com/nf-core/ampliseq/2.9.0/modules/local/dada2_taxonomy.nf -o modules/local/dada2/dada2_taxonomy.nf

# Fix PICRUST (If it's also showing that error)
curl -L https://raw.githubusercontent.com/nf-core/ampliseq/2.9.0/modules/local/picrust.nf -o modules/local/picrust.nf


nf-core pipelines schema build
````





Once have the csv from group A run:

````
nextflow run main.nf --input path/to/groupA_samplesheet.csv -profile docker
git add .
git commit -m "Group B: Integrated DADA2, Taxonomy, and PICRUSt logic with relative paths"
````

