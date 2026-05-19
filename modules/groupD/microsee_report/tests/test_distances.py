"""Tests for distance-matrix alignment."""

import pytest

from report_generator.charts.distances import align_distance_matrix
from report_generator.models import DistanceMatrixResult


def test_align_distance_matrix_reorders_samples() -> None:
    dm = DistanceMatrixResult(
        samples=["B", "A", "C"],
        matrix=[
            [0.0, 0.2, 0.3],
            [0.2, 0.0, 0.1],
            [0.3, 0.1, 0.0],
        ],
        n=3,
    )
    mat = align_distance_matrix(dm, ["A", "B", "C"])
    assert mat.shape == (3, 3)
    assert mat[0, 1] == pytest.approx(0.2)
    assert mat[1, 0] == pytest.approx(0.2)


def test_align_distance_matrix_missing_sample_raises() -> None:
    dm = DistanceMatrixResult(
        samples=["A", "B"],
        matrix=[[0.0, 0.5], [0.5, 0.0]],
        n=2,
    )
    with pytest.raises(ValueError, match="missing"):
        align_distance_matrix(dm, ["A", "X"])
