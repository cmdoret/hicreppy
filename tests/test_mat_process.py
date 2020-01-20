import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import scipy.stats as ss
import scipy.signal as sig
from scipy.stats import multivariate_normal
import hicreppy.utils.mat_process as hmp


### GENERATING SYNTHETIC DATA
def gauss_mat(meanx, meany, std, shape=(100, 100)):

    # create 2 kernels
    means = (meanx, meany)
    stds = np.eye(2) * std
    k = multivariate_normal(mean=means, cov=stds)

    # create a grid of (x,y) coordinates at which to evaluate the kernels
    x = np.linspace(-10, 10, shape[0])
    y = np.linspace(-10, 10, shape[1])
    xx, yy = np.meshgrid(x, y)

    # evaluate kernels at grid points
    xxyy = np.c_[xx.ravel(), yy.ravel()]
    zz = k.pdf(xxyy)

    # reshape and plot image
    g_mat = zz.reshape(shape)
    g_mat = sp.coo_matrix(g_mat)
    return g_mat

gauss1_mats = []
for mx in np.arange(-1.5, 1.5, 0.5):
    for my in np.arange(-1.5, 1.5, 0.5):
        # Only draw matrices with pattern in upper triangle
        if mx >= my:
            continue
        for sd in np.arange(0.3, 3, 0.3):
            m = gauss_mat(mx, my, sd)
            gauss1_mats.append(m)
GAUSS_MATS = ('signal', gauss1_mats)
N_CONTACTS = ('n_contacts', [2, 100, 10000])
PROP_CONTACTS = ('prop', [0, 0.1, 0.5, 0.8, 1])
MAT = sp.triu(sp.random(10000, 10000, density=0.05, format='coo'))

@pytest.mark.parametrize(*GAUSS_MATS)
def test_smooth(signal):
    """Check if smoothed matrix match output from scipy.signal"""
    h = 3
    km = 2 * h + 1
    kernel = np.ones((km, km)) / (km ** 2)
    corr_mat_sparse = hmp.smooth(signal, h).todense()
    # Use scipy result as base truth to compare smooth results
    corr_mat_scipy = np.zeros(signal.shape)
    corr_mat_scipy[h:-h, h:-h] = sig.correlate2d(
        signal.todense(), kernel, "valid"
    )
    assert np.allclose(corr_mat_sparse, corr_mat_scipy)


def test_subsample_contacts_exceed():
    """Oversampling should result in value errors"""
    with pytest.raises(ValueError):
        hmp.subsample_contacts(MAT, 10e10)

@pytest.mark.parametrize(*PROP_CONTACTS)
def test_subsample_contacts_prop(prop):
    """Test sampling proportions of contacts"""
    sampled = hmp.subsample_contacts(
        MAT.tocoo(), int(prop * MAT.data.sum())
    )
    assert np.isclose(sampled.data.sum(), MAT.data.sum() * prop, rtol=0.1)


@pytest.mark.parametrize(*N_CONTACTS)
def test_subsample_contacts_count(n_contacts):
    """Test sampling raw contact counts"""
    sampled = hmp.subsample_contacts(MAT.tocoo(), n_contacts)
    assert np.isclose(sampled.data.sum(), n_contacts, rtol=0.1)


@pytest.mark.filterwarnings("ignore::Warning")
@pytest.mark.parametrize(*GAUSS_MATS)
def test_diag_trim(signal):
    """Check if trimming diagonals preserves shape and sets diagonals to zero."""
    # Sparse version
    for d in range(signal.shape[0]):
        trimmed = hmp.diag_trim(signal.todia(), d)
        diag_sums = [
            trimmed.diagonal(d).sum() for d in range(trimmed.shape[0])
        ]
        assert trimmed.shape == signal.shape
        assert np.sum(diag_sums[d + 1 :]) == 0
    # Dense version
    signal = signal.toarray()
    for d in range(signal.shape[0]):
        trimmed = hmp.diag_trim(signal, d)
        diag_sums = [
            trimmed.diagonal(d).sum() for d in range(trimmed.shape[0])
        ]
        assert trimmed.shape == signal.shape
        assert np.sum(diag_sums[d + 1 :]) == 0


def test_set_mat_diag():
    """Test if dense matrix's nth diagonal is set to appropriate value"""
    dense_mat = np.zeros((10, 10))
    # The function operates in place
    hmp.set_mat_diag(dense_mat, diag=1, val=4)
    assert np.all(dense_mat.diagonal(1) == 4)
