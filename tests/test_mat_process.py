import numpy as np
import pytest
import scipy.sparse as sp
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
    points_x = np.linspace(-10, 10, shape[0])
    points_y = np.linspace(-10, 10, shape[1])
    grid_x, grid_y = np.meshgrid(points_x, points_y)

    # evaluate kernels at grid points
    grid = np.c_[grid_x.ravel(), grid_y.ravel()]
    gaussian = k.pdf(grid)

    # reshape and plot image
    g_mat = gaussian.reshape(shape)
    g_mat = sp.coo_matrix(g_mat)
    return g_mat


GAUSS1_MATS = []
for mx in np.arange(-1.5, 1.5, 0.5):
    for my in np.arange(-1.5, 1.5, 0.5):
        # Only draw matrices with pattern in upper triangle
        if mx >= my:
            continue
        for sd in np.arange(0.3, 3, 0.3):
            m = gauss_mat(mx, my, sd)
            GAUSS1_MATS.append(m)
GAUSS_MATS = ("signal", GAUSS1_MATS)
N_CONTACTS = ("n_contacts", [2, 100, 10000])
PROP_CONTACTS = ("prop", [0, 0.1, 0.5, 0.8, 1])
MAT = sp.triu(sp.random(10000, 10000, density=0.05, format="coo"))


@pytest.mark.parametrize(*GAUSS_MATS)
def test_smooth(signal):
    """Check if smoothed matrix match output from scipy.signal"""
    h_value = 3
    ker_m = 2 * h_value + 1
    kernel = np.ones((ker_m, ker_m)) / (ker_m**2)
    corr_mat_sparse = hmp.smooth(signal, h_value).todense()
    # Use scipy result as base truth to compare smooth results
    corr_mat_scipy = np.zeros(signal.shape)
    corr_mat_scipy[h_value:-h_value, h_value:-h_value] = sig.correlate2d(
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
    sampled = hmp.subsample_contacts(MAT.tocoo(), int(prop * MAT.data.sum()))
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
    for diag in range(signal.shape[0]):
        trimmed = hmp.diag_trim(signal.todia(), diag)
        diag_sums = [trimmed.diagonal(d).sum() for d in range(trimmed.shape[0])]
        assert trimmed.shape == signal.shape
        assert np.sum(diag_sums[diag + 1 :]) == 0
    # Dense version
    signal = signal.toarray()
    for d in range(signal.shape[0]):
        trimmed = hmp.diag_trim(signal, d)
        diag_sums = [trimmed.diagonal(d).sum() for d in range(trimmed.shape[0])]
        assert trimmed.shape == signal.shape
        assert np.sum(diag_sums[d + 1 :]) == 0


def test_set_mat_diag():
    """Test if dense matrix's nth diagonal is set to appropriate value"""
    dense_mat = np.zeros((10, 10))
    # The function operates in place
    hmp.set_mat_diag(dense_mat, diag=1, val=4)
    assert np.all(dense_mat.diagonal(1) == 4)
