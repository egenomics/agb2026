FROM python:3.11-slim

LABEL org.opencontainers.image.title="MicroSee Report Generator" \
      org.opencontainers.image.description="Python dependencies for the MicroSee Nextflow process. The report_generator/ code is NOT baked in — it is staged from the Nextflow work directory at runtime." \
      org.opencontainers.image.source="https://github.com/agb2026/agb2026"

RUN pip install --no-cache-dir \
    "pandas>=2.2,<3" \
    "numpy>=1.26,<3" \
    "pydantic>=2.7,<3"

WORKDIR /work
