"""
Cluster algorithms for phylogenetic approaches in linguistics.
"""
from itertools import combinations, product
from statistics import mean, median

from pylocluster.util import squareform, check_language_names


SCORES = {
    "upgma": mean,
    "average": mean,
    "single": min,
    "complete": max,
    "median": median,
}


def flat_linkage_recursive(clusters, matrix, threshold, linkage):
    if len(clusters) == 1:
        return
    scores, indices = [], []
    for (i, valA), (j, valB) in combinations(clusters.items(), r=2):
        score = []
        for vA, vB in product(valA, valB):
            score.append(matrix[vA][vB])
        scores.append(linkage(score))
        indices.append((i, j))
    minimum = min(scores)
    if minimum <= threshold:
        idxA, idxB = indices[scores.index(minimum)]
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return flat_linkage_recursive(clusters, matrix, threshold, linkage)


def flat_linkage(matrix, taxa=None, threshold=0.5, method="upgma", revert=True):
    clusters = dict([(i, [i]) for i in range(len(matrix))])
    out, tree = {}, []
    flat_linkage_recursive(clusters, matrix, threshold, SCORES[method])
    if taxa:
        for key in clusters:
            clusters[key] = [taxa[i] for i in clusters[key]]
    if revert:
        for idx, key in enumerate(clusters):
            for i in clusters[key]:
                out[i] = idx + 1
        return out
    return [set(cluster) for cluster in clusters.values()]


def linkage_recursive(clusters, matrix, tree, branches, linkage):
    """
    Recursive base version of linkage clustering.
    """
    if len(clusters) == 1:
        return

    scores, indices = [], []
    for (i, valA), (j, valB) in combinations(clusters.items(), r=2):
        score = []
        for vA, vB in product(valA, valB):
            score.append(matrix[vA][vB])
        scores.append(linkage(score))
        indices.append((i, j))
    minimum = min(scores)

    idxNew = max(clusters) + 1
    idxA, idxB = indices[scores.index(minimum)]

    bA = minimum / 2 - branches[idxA]
    bB = minimum / 2 - branches[idxB]

    branches[idxNew] = minimum / 2
    clusters[idxNew] = clusters[idxA] + clusters[idxB]

    del clusters[idxA]
    del clusters[idxB]

    tree.append([idxA, idxB, bA, bB])

    return linkage_recursive(clusters, matrix, tree, branches, linkage)


def linkage(matrix, taxa=None, method="upgma", distances=True):
    """
    Carry out a linkage clustering analysis.
    """
    formatter = "({0}:{2:.2f},{1}:{3:.2f})" if distances else "({0},{1})"
    taxa = check_language_names(taxa) or [
        "t_" + str(i + 1) for i in range(len(matrix[0]))
    ]
    clusters = dict([(i, [i]) for i in range(len(taxa))])
    branches = dict([(i, 0) for i in range(len(taxa))])
    newick = dict([(i, taxa[i]) for i in range(len(taxa))])
    tree = []

    linkage_recursive(clusters, matrix, tree, branches, SCORES[method])

    # create different output, depending on the options for the inclusion of
    # distances or topology only
    for i, (a, b, c, d) in enumerate(tree):
        newick[len(taxa) + i] = formatter.format(newick[a], newick[b], c, d)

    return newick[max(newick.keys())] + ";"


def neighbor_recursive(clusters, matrix, tree_matrix, constant_matrix, tracer):
    """
    Internal implementation of the neighbor-joining algorithm.
    """
    # terminate when the dictionary is of length 2
    if len(clusters) == 2:
        idxA, idxB = 0, 1
        idxNew = max(tracer.values()) + 1
        tracer[tuple(clusters[idxA] + clusters[idxB])] = idxNew
        sAX = matrix[idxA][idxB] / 2
        sBX = matrix[idxA][idxB] / 2
        tree_matrix.append(
            (tracer[tuple(clusters[idxA])], tracer[tuple(clusters[idxB])], sAX, sBX)
        )
        # join the clusters according to the index
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return

    # create the constant matrix when the process starts
    N = len(matrix)

    # determine the average scores (divergence r)
    averages = []
    for line in matrix:
        averages.append(sum(line) / (N - 2.0))

    # create the new matrix
    new_matrix = [[cell for cell in line] for line in matrix]
    for i, j in combinations(range(len(matrix)), r=2):
        new_score = matrix[i][j] - averages[i] - averages[j]
        new_matrix[i][j] = new_score
        new_matrix[j][i] = new_score

    # determine the minimal score
    scores, indices = [], []
    for i, j in combinations(clusters.keys(), r=2):
        scores.append(new_matrix[i][j])
        indices.append((i, j))

    minimum = min(scores)
    idxA, idxB = indices[scores.index(minimum)]

    # check for the average of the clusters
    vals = []
    for i, j in product(clusters[idxA], clusters[idxB]):
        vals.append(constant_matrix[i][j])
    tmp_score = sum(vals) / len(vals)

    # append the indices to the tree matrix
    sAX = matrix[idxA][idxB] / 2.0 + (averages[idxA] - averages[idxB]) / 2
    sBX = matrix[idxA][idxB] - sAX
    tree_matrix.append(
        (tracer[tuple(clusters[idxA])], tracer[tuple(clusters[idxB])], sAX, sBX)
    )

    # create the new index for the tracer
    idxNew = max(tracer.values()) + 1
    tracer[tuple(clusters[idxA] + clusters[idxB])] = idxNew

    clusters[idxA] += clusters[idxB]
    del clusters[idxB]

    new_clusters = {}
    for i, key in enumerate(sorted(clusters.keys())):
        new_clusters[i] = clusters[key]

    new_matrix = []
    dist_ab = matrix[idxA][idxB]
    for (i, a), (j, b) in combinations(enumerate(sorted(clusters.keys())), r=2):
        if i < j and a != idxA and b != idxA:
            new_matrix.append(matrix[a][b])
        elif i < j and a == idxA:
            dist_a = matrix[idxA][b]
            dist_b = matrix[idxB][b]
            new_matrix.append(((dist_a + dist_b) - dist_ab) / 2.0)
        elif i < j and b == idxA:
            dist_a = matrix[idxA][a]
            dist_b = matrix[idxB][a]
            new_matrix.append(((dist_a + dist_b) - dist_ab) / 2.0)

    # get values of new_clusters into clusters
    clusters = {}
    for key, val in new_clusters.items():
        clusters[key] = val

    # return score
    return neighbor_recursive(
        clusters, squareform(new_matrix), tree_matrix, constant_matrix, tracer
    )


def neighbor(matrix, taxa=None, distances=True):
    """
    Function clusters data according to the Neighbor-Joining algorithm \
    (:evobib:`Saitou1987`).
    """
    clusters = dict([(i, [i]) for i in range(len(taxa))])
    formatter = "({0}:{2:.2f},{1}:{3:.2f})" if distances else "({0},{1})"
    taxa = check_language_names(taxa) or [
        "t_" + str(i + 1) for i in range(len(matrix[0]))
    ]
    tracer = dict([(tuple([a]), b[0]) for (a, b) in clusters.items()])
    newick = dict([(i, taxa[i]) for i in range(len(matrix[0]))])

    tree = []

    neighbor_recursive(
        clusters,
        matrix,
        tree,
        [[c for c in l] for l in matrix],
        dict([(tuple([a]), b[0]) for (a, b) in clusters.items()]),
    )

    # create different output, depending on the options for the inclusion of
    # distances or topology only
    for i, (a, b, c, d) in enumerate(tree):
        newick[len(taxa) + i] = formatter.format(newick[a], newick[b], c, d)
    newick_string = newick[max(newick.keys())] + ";"
    return newick_string
