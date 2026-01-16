import numpy as np
from scipy import stats


def compare_independent_r2(r2_1, n1, r2_2, n2):
    print(f"Comparing R2 values: r2_1={r2_1}, n1={n1}, r2_2={r2_2}, n2={n2}")
    # Convert R-squared to correlation coefficient (r)
    r1 = np.sqrt(r2_1)
    r2 = np.sqrt(r2_2)

    # Fisher's r-to-z transformation
    z1 = np.arctanh(r1)
    z2 = np.arctanh(r2)

    # Standard error of the difference
    se_diff = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_observed = (z1 - z2) / se_diff

    # Determine the p-value (two-tailed)
    p_value = 2 * (1 - stats.norm.cdf(np.abs(z_observed)))
    print(f"Transformations: r1={r1:.4f}, r2={r2:.4f}, z1={z1:.4f}, z2={z2:.4f}, se_diff={se_diff:.4f}")
    print(f"Z observed: {z_observed:.4f}, p-value: {p_value:.4g}")
    return z_observed, p_value


n = 256
compare_independent_r2(0.95, n, 0.85, n)
compare_independent_r2(0.79, n, 0.96, n)
