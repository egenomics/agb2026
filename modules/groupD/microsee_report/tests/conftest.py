"""Shared pytest hooks for microsee_report."""

from __future__ import annotations

import sys


def pytest_configure(config) -> None:
    markexpr = config.getoption("markexpr", default="")
    if markexpr and "integration" in markexpr and "not integration" not in markexpr:
        sys.stderr.write(
            "\n[pytest] Running integration tests (CLI cohort HTML + in-process patient). "
            "Expect ~1–3 minutes, not 20.\n\n"
        )
