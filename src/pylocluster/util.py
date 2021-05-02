"""
Utility functions for the clustering package.
"""
from math import sqrt


def squareform(x):
    """
    A simplified version of the :py:func:`scipy.spatial.distance.squareform` \
    function.

    Parameters
    ----------

    x : list
        The one-dimensional flat representation of a symmetrix distance matrix.

    Returns
    -------
    matrix : list
        The two-dimensional redundant representation of a symmetric distance matrix.
    """
    l = len(x)

    # calculate the length of the square
    s = int(sqrt(2 * l) + 1)

    out = [[0.0 for i in range(s)] for j in range(s)]

    k = 0
    for i in range(s):
        for j in range(s):
            if i < j:
                out[i][j] = x[k]
                out[j][i] = x[k]
                k += 1
    return out


def check_language_names(taxa):
    for taxon in taxa or []:
        if [x for x in taxon if x in "():;, "]:
            raise ValueError(
                "No brackets, colons, spaces, and commas allowed for doculect names"
            )
    return taxa
