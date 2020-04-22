# hicreppy
**cmdoret**


[![PyPI version](https://badge.fury.io/py/hicreppy.svg)](https://badge.fury.io/py/hicreppy)
[![Build Status](https://travis-ci.com/cmdoret/hicreppy.svg)](https://travis-ci.com/cmdoret/hicreppy)
[![codecov](https://codecov.io/gh/cmdoret/hicreppy/branch/master/graph/badge.svg)](https://codecov.io/gh/cmdoret/hicreppy)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/cmdoret/hicreppy.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cmdoret/hicreppy/context:python)

This is a python reimplementation of hicrep's algorithm with added support for sparse matrices (in .cool format). 

hicrep measures similarity between Hi-C samples by computing a stratum-adjusted correlation coefficient (SCC). In this implementation, the SCC is computed separately for each chromosome and the chromosome length-weighted average of SCCs is computed.

hicrep is published at:
> HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng Yue, Qunhua Li, 2017, Genome Research, doi: 10.1101/gr.220640.117

The original implementation, in R can be found at https://github.com/MonkeyLB/hicrep

### Installation

You can install the package using pip:

```bash
pip install --user hicreppy
```

### Usage

To find the optimal value for smoothing parameter h, you can use the htrain subcommand:

```

Usage: hicreppy htrain [OPTIONS] COOL1 COOL2

  Find the optimal value for smoothing parameter h. The optimal h-value is
  printed to stdout. Run informations are printed to stderr.

Options:
  -r, --h-max INTEGER     Maximum value of the smoothing parameter h to
                          explore. All consecutive integer values from 0 to
                          this value will be tested.  [default: 10]
  -m, --max-dist INTEGER  Maximum distance at which to compute the SCC, in
                          basepairs.  [default: 100000]
  -b, --blacklist TEXT    Exclude those chromosomes in the analysis. List of
                          comma-separated chromosome names.
  -w, --whitelist TEXT    Only include those chromosomes in the analysis. List
                          of comma-separated chromosome names.
  --help                  Show this message and exit.

```

To compute the SCC between two matrices, use the scc subcommand. The optimal h value obtained with htrain should be provided to the flag `-v`:

```

Usage: hicreppy scc [OPTIONS] COOL1 COOL2

  Compute the stratum-adjusted correlation coefficient for input matrices

Options:
  -v, --h-value INTEGER    Value of the smoothing parameter h to use. Should
                           be an integer value >= 0.  [default: 10]
  -m, --max-dist INTEGER   Maximum distance at which to compute the SCC, in
                           basepairs.  [default: 100000]
  -s, --subsample INTEGER  Subsample contacts from both matrices to target
                           value. Leave to 0 to disable subsampling.
                           [default: 0]
  -b, --blacklist TEXT     Exclude those chromosomes in the analysis. List of
                           comma-separated chromosome names.
  -w, --whitelist TEXT     Only include those chromosomes in the analysis.
                           List of comma-separated chromosome names.
  --help                   Show this message and exit.
```

When running multiple pairwise comparisons, compute the optimal h value once between two highly similar samples and reuse the h value for all `scc` commands

### Contributing

All contributions are welcome. We use the [numpy standard](https://numpydoc.readthedocs.io/en/latest/format.html) for docstrings when documenting functions.

The code formatting standard we use is [black](https://github.com/psf/black), with --line-length=79 to follow PEP8 recommendations. We use `pytest` with the `pytest-doctest` and `pytest-pylint` plugins as our testing framework. Ideally, new functions should have associated unit tests, placed in the `tests` folder.

To test the code, you can run:

```bash
pytest --doctest-modules --pylint --pylint-error-types=EF --pylint-rcfile=.pylintrc hicreppy tests
```
