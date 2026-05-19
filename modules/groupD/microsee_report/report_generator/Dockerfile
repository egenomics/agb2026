FROM python:3.11-slim

LABEL org.opencontainers.image.title="MicroSee Report Generator" \
      org.opencontainers.image.description="Python dependencies for the MicroSee Nextflow process. \
The report_generator/ code is NOT baked in — it is staged from the Nextflow work directory at runtime." \
      org.opencontainers.image.source="https://github.com/agb2026/agb2026"

# Pin exact versions for reproducible builds.
# Regenerate: docker run --rm <this-image> pip freeze
RUN pip install --no-cache-dir \
    "pandas==2.3.3" \
    "numpy==2.4.1" \
    "pydantic==2.13.4"

# Smoke-test that all imports work before pushing the image.
RUN python -c "import pandas, numpy, pydantic; print('deps ok')"

WORKDIR /work
