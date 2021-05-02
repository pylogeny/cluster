from pylocluster.cluster import *
from pylocluster.util import squareform, check_language_names

def test_flat_linkage():
    matrix = squareform([0.5, 0.67, 0.8, 0.2, 0.4, 0.7, 0.6, 0.8, 0.8, 0.3])
    taxa = ['G', 'S', 'I', 'E', 'D']
    for score in SCORES:
        groups = flat_linkage(matrix, taxa=taxa, threshold=0.5, revert=True)
        assert groups['G'] == groups['D']

    matrix = squareform([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
    for score in SCORES:
        groups = flat_linkage(matrix, taxa=taxa, threshold=0.3, revert=False)
        assert len(groups) == 1


def test_linkage():
    matrix = squareform([0.5, 0.67, 0.8, 0.2, 0.4, 0.7, 0.6, 0.8, 0.8, 0.3])
    taxa = ['G', 'S', 'I', 'E', 'D']
    for score in SCORES:
        tree = linkage(matrix, taxa=taxa, method=score)
        tree2 = linkage(matrix, taxa=None, method=score, distances=False)
        assert 't_1' in tree2
        assert 'G' in tree


def test_neighbor():
    matrix = squareform([0.5, 0.67, 0.8, 0.2, 0.4, 0.7, 0.6, 0.8, 0.8, 0.3])
    taxa = ['G', 'S', 'I', 'E', 'D']
    tree = neighbor(matrix, taxa=taxa)

    matrix = squareform([0.5, 0.67, 0.8, 0.2, 0.4, 0.7, 0.6, 0.8, 0.8, 0.3])
    taxa = ['G', 'S', 'I', 'E', 'D']
    tree = neighbor(matrix, taxa=taxa)


def test_check_language_names():
    try:
        a = 1
        check_language_names(["bla", "bl,"])
        a = 2
    except ValueError:
        assert a == 1

