# pretsa_spatial.py
import numpy as np
from numba import njit
from patsy import dmatrix
from scipy.stats import f as fdist

def _bs(x, df):
    B = dmatrix(f"bs(x, df={df}, include_intercept=False) - 1",
                {"x": np.asarray(x, dtype=np.float64)},
                return_type="dataframe").to_numpy()
    return np.ascontiguousarray(B, dtype=np.float64)

def _tp_basis(row_vals, col_vals, df):
    X = _bs(row_vals, df)
    Y = _bs(col_vals, df)
    X = np.concatenate([np.ones((X.shape[0], 1)), X], axis=1)
    Y = np.concatenate([np.ones((Y.shape[0], 1)), Y], axis=1)
    cols = [X[:, i] * Y[:, j] for i in range(X.shape[1]) for j in range(Y.shape[1])]
    B = np.stack(cols, axis=1)
    keep = B.std(axis=0, ddof=1) > 0
    B = B[:, keep]
    B = np.concatenate([np.ones((B.shape[0], 1)), B], axis=1)
    return np.ascontiguousarray(B, dtype=np.float64)

@njit(cache=True)
def _pred(expr_GN, B_NP):
    tBB = B_NP.T @ B_NP
    inv_tBB = np.linalg.inv(tBB)
    return (expr_GN @ B_NP) @ inv_tBB @ B_NP.T

def spatialTest(expr, coord, knot=0, maxknotallowed=5):
    expr = np.asarray(expr, dtype=np.float64)
    row_vals = np.asarray(coord["row"], dtype=np.float64)
    col_vals = np.asarray(coord["col"], dtype=np.float64)
    df = int(knot) + 3
    B = _tp_basis(row_vals, col_vals, df)
    pred = _pred(expr, B)

    SSE = np.sum((expr - pred) ** 2, axis=1)
    mean_row = np.mean(expr, axis=1, keepdims=True)
    SST = np.sum((expr - mean_row) ** 2, axis=1)

    fstat = ((SST - SSE) / (B.shape[1] - 1)) / (SSE / (B.shape[0] - B.shape[1]))
    zero_mask = (np.sum(expr, axis=1) == 0)
    fstat[zero_mask] = 0.0

    df1 = B.shape[1] - 1
    df2 = B.shape[0] - B.shape[1]
    pval = fdist.sf(fstat, df1, df2)
    logpval = fdist.logsf(fstat, df1, df2)

    out = np.column_stack([logpval, pval, fstat])
    return out  # columns: logpval, pval, fstat

def spatialFit(expr, coord, knot=0, maxknotallowed=5):
    expr = np.asarray(expr, dtype=np.float64)
    row_vals = np.asarray(coord["row"], dtype=np.float64)
    col_vals = np.asarray(coord["col"], dtype=np.float64)
    df = int(knot) + 3
    B = _tp_basis(row_vals, col_vals, df)
    pred = _pred(expr, B)
    return pred
