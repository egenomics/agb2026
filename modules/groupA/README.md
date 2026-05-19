# GROUP A: Data Handling & Preprocessing

This folder contains all modules developed by Group A.

## Group Responsibilities
- Data input and validation
- Quality control
- Preprocessing and normalization
- Metadata integration

## Workflow

```mermaid
---
title: Workflow Group A
---
flowchart 

%% Sofware / tools = rect
%% Nextflow modules = st-rect
%% Database = cyl
%% Input = sl-rect
%% Output = lean
%% documentation = doc
%% decision = diamond
%% comments = brace
%% Manual input = lean-r


metadata@{shape: sl-rect, label: "Metadata"}
comment@{shape : braces, label: "Format metadata to match the required </br>format for the workflow" }
format@{shape: lin-doc, label: "Metadata Format"}

%%doc@{shape : doc, label: "Documentation"}
%%O@{ shape: lean-l, label: "Preprocessed Data" }
%%I@{ shape: sl-rect, label: "Raw Data" }

A@{ shape: circle, label: "Start" }
type@{ shape: diamond, label: "Data type" }
Qiita@{ shape: rect, label: "Qiita" }
prefetch@{ shape: trap-t, label: "Prefetch" }
NCBI@{ shape: cyl, label: "NCBI" }
fastq@{ shape: sl-rect, label: "SRA files" }
split@{ shape: st-rect, label: "FASTQ dump (PE)"}
split2@{ shape: st-rect, label: "FASTQ dump (SE)"}
se@{ shape: lean-r, label: "SE.fastq"}
qreport@{shape: lean-l, label: "FastQC report"}

pe_1@{shape: lean-l, label: "PE_1.fastq"}
pe_2@{shape: lean-r, label: "PE_2.fastq"}
FASTQC@{shape: st-rect, label: "FASTQC"}
CUTADAPT@{shape: st-rect, label: "Cutadapt"}
FASTQC2@{shape: st-rect, label: "FASTQC"}
MultiQC@{shape: st-rect, label: "MultiQC"}
trimmed@{shape: lean-l, label: "FASTQ"}


%% Diagram

A --> Qiita --> metadata


subgraph NC [ Retrieve Data ]
    direction LR
    metadata --> prefetch --> NCBI
end

NCBI --> fastq
comment --> format

subgraph Nextflow ["Nextflow Pipeline"]
    direction TB

    fastq & format --> type
    type -->|PE| split --> pe_1 & pe_2 --> FASTQC
    type -->|SE| split2 --> se --> FASTQC
    FASTQC --> CUTADAPT --> FASTQC2
    FASTQC2 --> MultiQC





end

split --- comment1@{shape: braces, label: "Split the strands?? <br> Alberto suggests" }
FASTQC ---> qreport


```

## Modules
<!-- Each group member should add their modules here -->

### Module Template
```
module_name/
├── main.nf          # Process definition
├── meta.yml         # Module metadata
└── tests/
    ├── main.nf.test      # Test definition
    └── main.nf.test.snap # Test snapshots
```

## Integration Points
This group outputs data consumed by:
- **Group B** - Receives preprocessed data

## Communication
- Document your outputs clearly
- Update the main workflow integration points
- Test your modules before integration

