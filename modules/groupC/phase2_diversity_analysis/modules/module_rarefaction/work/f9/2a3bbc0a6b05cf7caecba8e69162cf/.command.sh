#!/bin/bash -ue
cat > rarefaction_skipped_report.txt << REPORT
========================================================
RAREFACTION SKIPPED
========================================================
Date            : $(date '+%Y-%m-%d %H:%M:%S')
Sample count    : 75
Minimum required: 100

Automatic rarefaction threshold selection was not
performed because the number of samples (75)
did not reach the minimum required (100).

Action required: accumulate more samples and re-run
the pipeline once the threshold is met.
========================================================
REPORT
