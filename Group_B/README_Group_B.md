# Group_B folder's organization

````
alphaflow_backup/            <-- If we need some files or folders to move from here to inside Group_B folder
Group_B/
├── biodb/                     <-- This folder won't be found on GitHub due to have huge weight (it will be on our shared drive :)
├── modules
├── make_samplesheet.sh        <-- The new automation script
├── main.nf                    <-- Entry point
├── nextflow.config            <-- Configured for cluster (slurm, etc.)
├── subworkflows/              <-- The 3 main folders
│   ├── qc_checks/             <-- (New: Created by you)
│   ├── taxonomic_profiling/   <-- (New: Created by you)
│   ├── functional_annotation/ <-- (New: Created by you)
│   ├── local/utils_.../main.nf          <-- (Keep: The template's helpers)
│   └── nf-core/utils_.../main.nf        <-- (Keep: The template's helpers)
└── assets/
    └── samplesheet.csv        <-- Automatically generated here
````
# About each subworkflow

## `qc_checks`/ (Quality Control)
This folder handles the "health check" of the data. It’s where you ensure Group A’s cleaning actually worked.

What’s inside:

- `fastqc.nf`: A module to check read quality.

- `multiqc.nf`: A module that compiles all logs into one final HTML report.

- The Output: The "Quality Report" that tells you if the data is good enough to trust.

## `taxonomic_profiling`/ (The "Who is there?")
It contains the logic for identifying the microbes.

What’s inside:

- `dada2_denoising.nf`: The logic for turning raw reads into ASVs.

- `dada2_taxonomy.nf`: The logic for comparing ASVs against the SILVA database.

- `kraken2.nf`: The ultra-fast k-mer classifier for a "second opinion."

- The Output: The ASV Table and the Taxonomy Table.

## `functional_annotation`/ (The "What are they doing?")

This folder takes the results from the profiling and predicts the biological impact.

What’s inside:

- `picrust2.nf`: The module that predicts metabolic pathways.

- `pathway_analysis.nf`: (Optional) Scripts to turn those predictions into readable charts (like MetaCyc pathways).

- The Output: A table of Metabolic Pathways (e.g., how much Insulin or Butyrate the community can produce).

## `main.nf` --> master one :)
One master `main.nf` at the top that "calls" the code inside these folders.

# How it fits the modular folders

Now that  have the three main folders (qc_checks, taxonomic_profiling, functional_annotation), here is how the data flows:

- The Bash Script: Creates assets/samplesheet.csv.

- Top-level main.nf: Reads that CSV and creates a "Channel."

- he Folders: * The qc_checks folder gets the raw paths first.

- The output of taxonomic_profiling (the ASVs) becomes the input for functional_annotation.

# TO DO LIST:

1. The Bash Script Logic (`make_samplesheet.sh`)
Your script should live in your alphaflow/bin/ folder or the root. It needs to find the files and handle the "Forward" (R1) and "Reverse" (R2) pairing. BE CAREFUL BECAUSE THE GROUP A SAID THAT THEY PERHAPS WILL BE WORKING ON SINGLE END DATA!! SO TAKE INTO ACCOUNT THIS!!!!!!

What the script should do:

- Search: Use find or ls in the cluster's data directory.

- Match: Identify pairs based on the _R1_ and _R2_ naming convention.

- Format: Output the results into a Comma-Separated format.

- Permissions: Make sure the script has execution rights (`chmod +x make_samplesheet.sh`) and that the user running Nextflow has "Read" access to the folder where Group A stored the files.