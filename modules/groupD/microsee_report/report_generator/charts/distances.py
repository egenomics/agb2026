"""charts/distances.py — distance matrices, PCoA ordination, and hierarchical clustering."""

from __future__ import annotations
import math
import numpy as np


def rows_to_ab(rows: list[dict], taxa: list[str]) -> np.ndarray:
    return np.array([[float(r.get(t) or 0) for t in taxa] for r in rows], dtype=float)


def bray_curtis_matrix(ab: np.ndarray) -> np.ndarray:
    n = ab.shape[0]
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.sum(np.abs(ab[i] - ab[j]))
            den = np.sum(ab[i] + ab[j])
            d = float(num / den) if den > 0 else 0.0
            mat[i, j] = mat[j, i] = d
    return mat


def jaccard_matrix(ab: np.ndarray, threshold: float = 5.0) -> np.ndarray:
    n = ab.shape[0]
    totals = ab.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    pres = (ab / totals * 100) >= threshold
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            either = int(np.sum(pres[i] | pres[j]))
            both   = int(np.sum(pres[i] & pres[j]))
            d = 1.0 - both / either if either > 0 else 0.0
            mat[i, j] = mat[j, i] = d
    return mat


def pcoa(dist_mat: np.ndarray) -> tuple[list[float], list[float], float, float]:
    n = dist_mat.shape[0]
    if n < 3:
        return [0.0] * n, [0.0] * n, 0.0, 0.0
    D2 = dist_mat ** 2
    row_mean   = D2.mean(axis=1, keepdims=True)
    col_mean   = D2.mean(axis=0, keepdims=True)
    grand_mean = D2.mean()
    B = -0.5 * (D2 - row_mean - col_mean + grand_mean)
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    pos   = eigenvalues > 0
    total = float(eigenvalues[pos].sum()) if pos.any() else 1.0
    lam1  = max(float(eigenvalues[0]), 0.0)
    lam2  = max(float(eigenvalues[1]), 0.0) if len(eigenvalues) > 1 else 0.0
    cx    = (eigenvectors[:, 0] * math.sqrt(lam1)).tolist()
    cy    = (eigenvectors[:, 1] * math.sqrt(lam2)).tolist() if len(eigenvalues) > 1 else [0.0] * n
    pct1  = round(lam1 / total * 100, 1) if total else 0.0
    pct2  = round(lam2 / total * 100, 1) if total else 0.0
    return cx, cy, pct1, pct2


def average_linkage_dendrogram(mat: np.ndarray) -> list:
    """Average-linkage clustering. Returns list of (xs, ys) segment tuples."""
    n = mat.shape[0]
    D: dict[int, dict[int, float]] = {
        i: {j: float(mat[i, j]) for j in range(n) if j != i}
        for i in range(n)
    }
    pos    = {i: float(i) for i in range(n)}
    height = {i: 0.0 for i in range(n)}
    size   = {i: 1 for i in range(n)}
    active = set(range(n))
    next_id = n
    segments: list = []

    for _ in range(n - 1):
        arr = sorted(active)
        min_d = float("inf"); ci = cj = -1
        for a in range(len(arr)):
            for b in range(a + 1, len(arr)):
                d = D.get(arr[a], {}).get(arr[b], float("inf"))
                if d < min_d:
                    min_d, ci, cj = d, arr[a], arr[b]
        if ci < 0:
            break

        pi, pj = pos[ci], pos[cj]
        hi, hj = height[ci], height[cj]
        ni, nj = size[ci], size[cj]
        segments.append(([hi, min_d, min_d, hj], [pi, pi, pj, pj]))

        D[next_id] = {}
        for k in active:
            if k in (ci, cj):
                continue
            dik = D.get(ci, {}).get(k, D.get(k, {}).get(ci, 0.0))
            djk = D.get(cj, {}).get(k, D.get(k, {}).get(cj, 0.0))
            d   = (dik * ni + djk * nj) / (ni + nj)
            D[next_id][k] = d
            D.setdefault(k, {})[next_id] = d

        pos[next_id]    = (pi * ni + pj * nj) / (ni + nj)
        height[next_id] = min_d
        size[next_id]   = ni + nj
        active.discard(ci); active.discard(cj); active.add(next_id)
        next_id += 1

    return segments
