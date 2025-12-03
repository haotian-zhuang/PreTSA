# pretsa_temporal.py
import numpy as np
from patsy import dmatrix
from scipy.stats import f as fdist, gamma as gamma_dist

def _bs(x, df):
    B = dmatrix(f"bs(x, df={df}, include_intercept=False) - 1",
                {"x": np.asarray(x, dtype=np.float64)},
                return_type="dataframe").to_numpy()
    return np.ascontiguousarray(B, dtype=np.float64)

def _design(pseudotime, df):
    B = _bs(pseudotime, df)
    keep = B.std(axis=0, ddof=1) > 0
    B = B[:, keep]
    B = np.concatenate([np.ones((B.shape[0], 1)), B], axis=1)
    return np.ascontiguousarray(B, dtype=np.float64)

def _pred(expr_GN, B_NP):
    tBB = B_NP.T @ B_NP
    coef = np.linalg.solve(tBB, B_NP.T @ expr_GN.T)
    return (B_NP @ coef).T

def Calfstat(expr, pseudotime, knot=0, maxknotallowed=10):
    expr = np.asarray(expr, dtype=np.float64)
    pt = np.asarray(pseudotime, dtype=np.float64)
    df = int(knot) + 3
    B = _design(pt, df)
    pred = _pred(expr, B)
    SSE = np.sum((expr - pred) ** 2, axis=1)
    mean_row = np.mean(expr, axis=1, keepdims=True)
    SST = np.sum((expr - mean_row) ** 2, axis=1)
    fstat = ((SST - SSE) / (B.shape[1] - 1)) / (SSE / (B.shape[0] - B.shape[1]))
    zero_mask = (np.sum(expr, axis=1) == 0)
    fstat[zero_mask] = 0.0
    return fstat

def temporalTestFixed(expr, pseudotime, knot=0, maxknotallowed=10):
    expr = np.asarray(expr, dtype=np.float64)
    pt = np.asarray(pseudotime, dtype=np.float64)
    df = int(knot) + 3
    B = _design(pt, df)
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
    out = np.column_stack([fstat, pval, logpval])
    return out  # columns: fstat, pval, logpval

def temporalTest(expr, pseudotime, pseudotime_permute=None, knot=0, maxknotallowed=10):
    if pseudotime_permute is None:
        return temporalTestFixed(expr=expr, pseudotime=pseudotime, knot=knot, maxknotallowed=maxknotallowed)
    else:
        fstat_ori = Calfstat(expr=expr, pseudotime=pseudotime, knot=knot, maxknotallowed=maxknotallowed)
        perms = []
        for pt in pseudotime_permute:
            perms.append(Calfstat(expr=expr, pseudotime=np.asarray(pt), knot=knot, maxknotallowed=maxknotallowed))
        fstat_perm = np.column_stack(perms) if len(perms) > 0 else np.empty((fstat_ori.shape[0], 0))
        pval_emp = (np.sum(fstat_perm >= fstat_ori[:, None], axis=1) + 1.0) / (fstat_perm.shape[1] + 1.0)
        pval_par = np.full_like(pval_emp, np.nan, dtype=np.float64)
        for i in range(fstat_perm.shape[0]):
            x = fstat_perm[i, :]
            if x.size >= 2 and np.all(np.isfinite(x)) and np.nanmax(x) > 0:
                try:
                    k, loc, theta = gamma_dist.fit(x, floc=0.0)
                    pval_par[i] = gamma_dist.sf(fstat_ori[i], a=k, loc=loc, scale=theta)
                except Exception:
                    pval_par[i] = np.nan
        out = np.column_stack([pval_par, pval_emp, fstat_ori])
        return out  # columns: pval.parametric, pval.empirical, fstat.ori

def temporalFit(expr, pseudotime, knot=0, maxknotallowed=10):
    expr = np.asarray(expr, dtype=np.float64)
    pt = np.asarray(pseudotime, dtype=np.float64)
    df = int(knot) + 3
    B = _design(pt, df)
    pred = _pred(expr, B)
    return pred
