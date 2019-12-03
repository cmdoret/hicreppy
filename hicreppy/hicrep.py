import sys
import os
import hicreppy.utils.mat_process as cu
import numpy as np
import scipy.stats as ss
from scipy.sparse import SparseEfficiencyWarning
import warnings


def h_train(mat1, mat2, max_dist, h_max, whitelist=None, blacklist=None):
    """
    Find the optimal value for the smoothing parameter h.
    For each value of h in h_range, in ascending order, separate intrachromosomal
    sub-matrices of each chromosome. Each sub-matrix is subsampled to 10% contacts
    and the stratum adjusted correlation coefficient (SCC) is computed 10 times between
    the corresponding chromosome of both samples. The mean of those 10 repetitions is computed
    for each chromosome and their mean weighted by chromosome length is used. The optimal h
    value is defined as the smallest value for which the increment of scc is less than 0.01.


    Parameters
    ----------
    mat1 : cooler.Cooler
        First matrix to compare. A Cooler object returned when loading a cool file.
    mat2 : cooler.Cooler
        Second matrix to compare. A Cooler object returned when loading a cool file.
    max_dist : int
        Maximum distance at which ton consider interactions, in basepairs.
    h_max : int 
        The maximum value of the smoothing parameter h (neighbourhood size) to test.
        All values between 0 and h_max will be explored.
    whitelist : None or list of strs
        If given, only compare those chromosomes.
    blacklist : None or list of strs
        If given, do not compare those chromosomes.

    Returns
    -------
    int : 
        Optimal value for the smoothing parameter h.
    """
    if mat1.binsize != mat2.binsize:
        raise ValueError("Both matrices must be binned at the same resolution")

    max_bins = max_dist // mat1.binsize
    # Define chromosomes to scan
    # NOTE: chromosomes smaller than the kernel used for smoothing or the max_dist must
    # be removed.
    min_size = max((2 * h_max + 1) * mat1.binsize, max_dist)
    chromlist, chrom_lengths = make_chromlist(
        mat1, whitelist, blacklist, min_size=min_size
    )

    prev_scc = -np.inf
    for h, h_value in enumerate(range(h_max)):
        # Compute SCC values separately for each chromosome
        chrom_scc = np.zeros(len(chromlist))
        for c, chrom in enumerate(chromlist):
            sample_scc = np.zeros(10)
            chrom_1 = mat1.matrix(sparse=True, balance=False).fetch(chrom)
            chrom_2 = mat2.matrix(sparse=True, balance=False).fetch(chrom)
            # Sample 10% contacts and smooth 10 times for this chromosome
            for sample in range(10):
                sub_1 = cu.subsample_contacts(chrom_1, 0.1)
                sub_2 = cu.subsample_contacts(chrom_2, 0.1)
                smooth_1 = cu.smooth(sub_1, h_value)
                smooth_2 = cu.smooth(sub_2, h_value)
                sample_scc[sample] = get_scc(
                    smooth_1, smooth_2, max_bins=max_bins
                )
            # Use average SCC from 10 subsamples
            chrom_scc[c] = np.mean(sample_scc)
        print(
            f"Found SCC of {round(chrom_scc[c], 3)} with h={h_value}.",
            file=sys.stderr,
        )
        # Compute the genome SCC for this value of h using the weighted averge
        # of chromosomes SCC by their lengths. NaN values of SCC are not considered
        # This happens when comparing empty diagonals
        nan_scc_mask = ~np.isnan(chrom_scc)
        trunc_scc = chrom_scc[nan_scc_mask]
        trunc_lengths = np.array(chrom_lengths)[nan_scc_mask]
        curr_scc = np.average(trunc_scc, weights=trunc_lengths)
        # Check if SCC improvement is less than threshold
        if curr_scc - prev_scc < 0.01:
            break
        prev_scc = curr_scc
    if h_value == h_max:
        print(
            "Note: It's likely that your searching range is too"
            "narrow. Try to expand the range and rerun it.\n",
            file=sys.stderr,
        )
    # Return last h value with improvement >= 0.01
    return max(h - 1, 0)


def genome_scc(
    mat1, mat2, max_dist, h, subsample=0, whitelist=None, blacklist=None
):
    """Compute the Stratum-adjusted correlation coefficient (SCC) for the whole genome
    from cool files.

    Separate intrachromosomal sub-matrices of each chromosome. Compute the stratum adjusted
    correlation coefficient (SCC) between the corresponding chromosome of both samples.
    The mean of SCCs weighted by chromosome length is used. 

    Parameters
    ----------
    mat1 : cooler.Cooler
        First matrix to compare. A Cooler object returned when loading a cool file.
    mat2 : cooler.Cooler
        Second matrix to compare. A Cooler object returned when loading a cool file.
    max_dist : int
        Maximum distance at which ton consider interactions, in basepairs.
    h : int
        The smoothing parameter. A mean filter is used to smooth matrices, this value
        is the size of the filter.
    subsample : int
        The number of contacts to which matrices should be subsampled. When you plan to
        compare multiple matrices, it can be useful to subsample all of them to the
        same value to remove potential biases caused by different coverages. Set to 0 to
        disable subsampling.
    whitelist : None or list of strs
        If given, only compare those chromosomes.
    blacklist : None or list of strs
        If given, do not compare those chromosomes.

    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    if mat1.binsize != mat2.binsize:
        raise ValueError("Both matrices must be binned at the same resolution")

    max_bins = max_dist // mat1.binsize

    # Define chromosomes to scan
    # NOTE: chromosomes smaller than the kernel used for smoothing or the max_dist must
    # be removed.
    min_size = max((2 * h + 1) * mat1.binsize, max_dist)
    chromlist, chrom_lengths = make_chromlist(
        mat1, whitelist, blacklist, min_size=min_size
    )

    # Compute SCC values separately for each chromosome
    chrom_scc = np.zeros(len(chromlist))
    for c, chrom in enumerate(chromlist):
        chrom_1 = mat1.matrix(sparse=True, balance=False).fetch(chrom)
        chrom_2 = mat2.matrix(sparse=True, balance=False).fetch(chrom)
        # Sample 10% contacts and smooth 10 times for this chromosome
        if subsample > 0:
            chrom_1 = cu.subsample_contacts(chrom_1, subsample)
            chrom_2 = cu.subsample_contacts(chrom_2, subsample)
        smooth_1 = cu.smooth(chrom_1, h)
        smooth_2 = cu.smooth(chrom_2, h)

        chrom_scc[c] = get_scc(smooth_1, smooth_2, max_bins=max_bins)
    # Compute the genome SCC using the weighted averge of chromosomes
    # SCC by their lengths. NaN values of SCC are not considered
    # This happens when comparing empty diagonals
    nan_scc_mask = ~np.isnan(chrom_scc)
    chrom_scc = chrom_scc[nan_scc_mask]
    chrom_lengths = np.array(chrom_lengths)[nan_scc_mask]
    scc = np.average(chrom_scc, weights=chrom_lengths)
    return scc


def get_scc(mat1, mat2, max_bins):
    """
    Compute the stratum-adjusted correlation coefficient (SCC) between two
    Hi-C matrices up to max_dist. A Pearson correlation coefficient is computed
    for each diagonal in the range of 0 to max_dist and a weighted sum of those
    coefficients is returned.

    Parameters
    ----------
    mat1 : scipy.sparse.csr_matrix
        First matrix to compare.
    mat2 : scipy.sparse.csr_matrix
        Second matrix to compare.
    max_bins : int
        Maximum distance at which to consider, in bins.

    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
        trimmed_1 = cu.diag_trim(mat1.todia(), max_bins).tocsr()
        trimmed_2 = cu.diag_trim(mat2.todia(), max_bins).tocsr()
    corr_diag = np.zeros(len(range(max_bins)))
    weight_diag = corr_diag.copy()
    for d in range(max_bins):
        d1 = trimmed_1.diagonal(d)
        d2 = trimmed_2.diagonal(d)
        # Silence NaN warnings: this happens for empty diagonals and will
        # not be used in the end.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ss.PearsonRConstantInputWarning
            )
            # Compute raw pearson coeff for this diag
            corr_diag[d] = ss.pearsonr(d1, d2)[0]
        # Compute weight for this diag
        r2k = cu.vstrans(d1, d2)
        weight_diag[d] = len(d1) * r2k
    # Normalize weights
    weight_diag /= sum(weight_diag)

    # Weighted sum of coefficients to get SCCs
    scc = np.sum(corr_diag * weight_diag)

    return scc


def make_chromlist(c, whitelist, blacklist, min_size=None):
    """Given a cool object, a blacklist and whitelist of chromosomes, return the list
    of chromosomes to include in the analysis.

    Parameters
    ----------
    c : cooler.Cooler
        First matrix to compare. A Cooler object returned when loading a cool file.
    whitelist : None or list of strs
        If given, only compare those chromosomes.
    blacklist : None or list of strs
        If given, do not compare those chromosomes.
    min_size : int
        Chromosomes smaller than this value will be removed.

    Returns
    -------
    chromlist : list of strs
        Names of chromosomes to process
    chromlist : list of strs
        Lengths of chromosomes to process
    """
    chroms = c.chroms()[:].set_index("name")
    chromlist = chroms.index.values.tolist()
    if whitelist is not None:
        chromlist = whitelist
    if blacklist is not None:
        for black in blacklist:
            chromlist.remove(black)
    # Remove small chromosomes
    orig_chromlist = np.array(chromlist).copy()
    if min_size is not None:
        for chrom in orig_chromlist:
            if chroms.loc[chrom].length < min_size:
                chromlist.remove(chrom)
    n_removed = len(orig_chromlist) - len(chromlist)
    print(
        f"Removed {n_removed} chromosomes shorter than max_dist",
        file=sys.stderr,
    )
    chrom_lengths = chroms.loc[chromlist].length.values.tolist()
    return chromlist, chrom_lengths
