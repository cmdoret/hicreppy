import numpy as np
import scipy.sparse as sp
import hicreppy.hicrep as hcr

MAT1 = sp.triu(sp.random(10000, 10000, density=0.05, format="coo"))
MAT2 = sp.triu(sp.random(10000, 10000, density=0.05, format="coo"))


def test_get_scc():
    max_dist = 1000
    assert np.isclose(hcr.get_scc(MAT1, MAT1, max_bins=max_dist), 1, rtol=10e-3)
    assert hcr.get_scc(MAT1, sp.coo_matrix(MAT1.shape), max_bins=max_dist) == 0
    assert abs(hcr.get_scc(MAT1, MAT2, max_bins=max_dist)) <= 1
