import matplotlib.pyplot as plt
import os
import numpy as np
from .cfg import * 
# import subprocess
# subprocess.call(["git","clone","https://github.com/MartijnGosgens/validation_indices"])
from BregmanInitializer.init_cluster import frommembershipMatriceToVector, fromVectorToMembershipMatrice

## taken from sklearn
def gen_even_slices(n, n_packs, *, n_samples=None):
    """Generator to create `n_packs` evenly spaced slices going up to `n`.

    If `n_packs` does not divide `n`, except for the first `n % n_packs`
    slices, remaining slices may contain fewer elements.

    Parameters
    ----------
    n : int
        Size of the sequence.
    n_packs : int
        Number of slices to generate.
    n_samples : int, default=None
        Number of samples. Pass `n_samples` when the slices are to be used for
        sparse matrix indexing; slicing off-the-end raises an exception, while
        it works for NumPy arrays.

    Yields
    ------
    `slice` representing a set of indices from 0 to n.

    See Also
    --------
    gen_batches: Generator to create slices containing batch_size elements
        from 0 to n.

    Examples
    --------
    >>> from sklearn.utils import gen_even_slices
    >>> list(gen_even_slices(10, 1))
    [slice(0, 10, None)]
    >>> list(gen_even_slices(10, 10))
    [slice(0, 1, None), slice(1, 2, None), ..., slice(9, 10, None)]
    >>> list(gen_even_slices(10, 5))
    [slice(0, 2, None), slice(2, 4, None), ..., slice(8, 10, None)]
    >>> list(gen_even_slices(10, 3))
    [slice(0, 4, None), slice(4, 7, None), slice(7, 10, None)]
    """
    start = 0
    if n_packs < 1:
        raise ValueError("gen_even_slices got n_packs=%s, must be >=1" % n_packs)
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            if n_samples is not None:
                end = min(n_samples, end)
            yield range(start, end)
            start = end
