# PyloCluster 

[![codecov](https://codecov.io/gh/pylogeny/cluster/branch/master/graph/badge.svg)](https://codecov.io/gh/pylogeny/cluster)
[![PyPI](https://img.shields.io/pypi/v/pylocluster.svg)](https://pypi.org/project/pylocluster)

PyloCluster provides basic functionalities for distance-based clustering procedures in Python, including implementations of the Neighbor-joining and the UPGMA algorithm for phylogenetic reconstruction.

## Installation

```shell
$ pip install pylocluster
```

## Usage

The following examples requires the [python-newick](https://github.com/glottobank/python-newick) package.

```python
>>> from pylocluster import *
>>> from newick import loads
>>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])
>>> nwk = linkage(matrix, taxa=['G', 'S', 'I', 'E', 'D'], method='upgma')
>>> print(loads(nwk).ascii_art())
        ┌─S
    ┌───┤
    │   └─I
────┤
    │   ┌─E
    └───┤
        │   ┌─G
        └───┤
            └─D
```
