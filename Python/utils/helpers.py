import numpy as np


def is_base_pair(bs1: str, bs2: str, has_wobble: bool, has_uu: bool) -> bool:
    """
    Return whether two bases form a pair.

    Args:
        bs1 (str): the first base.
        bs2 (str): the second base.
        has_wobble (bool): whether G-U pair is allowed.
        has_uu (bool): whether U-U pair is allowed.

    Returns:
        bool: whether bs1 and bs2 forms a pair.
    """
    base_map = {
        'A': ['T', 'U'],
        'U': ['A'],
        'C': ['G'],
        'G': ['C'],
        'T': ['A']
    }
    if has_wobble:
        base_map['G'].append('U')
        base_map['U'].append('G')
    if has_uu:
        base_map['U'].append('U')

    if bs2.upper() in base_map[bs1.upper()]:
        return True
    return False


def filter_based_on_stem_loop_score(vertex_set: list, sl_min: float,
                                    sl_max: float) -> list:
    """
    Filter a given vertex set based on given parameters.

    Args:
        vertex_set (list): a list of vertex.
        sl_min (float): minimum stem-loop-score allowed in the filtered set.
        sl_max (float): maximum stem-loop-score allowed in the filtered set.

    Returns:
        list: a list of filtered vertex.
    """
    filtered_vertex_set = []
    count = 0
    for i in range(len(vertex_set)):
        v = vertex_set[i]
        if v.stem_loop_score <= sl_max and v.stem_loop_score >= sl_min:
            v.set_index(count)
            count += 1
            filtered_vertex_set.append(vertex_set[i])
    return filtered_vertex_set


def filter_based_on_distance_and_open_ended_modification(
        vertex_set: list, sl_min: float, sl_max: float, dist_min: int,
        dist_max: int, l: int) -> list:
    """
    Filter a given vertex set based on given parameters.

    Args:
        vertex_set (list): a list of vertex.
        sl_min (float): minimum stem-loop-score.
        sl_max (float): maximum stem-loop-score.
        dist_min (int): minimum distance of a stem.
        dist_max (int): maximum distance of a stem.
        l (int): length of the given rna sequence.

    Returns:
        list: a list of filtered possible vertex.
    """
    filtered_vertex_set = []
    ct = 0
    for i in range(len(vertex_set)):
        v = vertex_set[i]

        if v.distance > l / 2:
            new_distance = -v.distance + l + 2 * v.stem_length - 2
            new_stem_loop_score = new_distance / v.stem_length
            if new_stem_loop_score < 3:
                v.distance = new_distance
                v.set_stem_loop_score()
                v.set_index(ct)
                ct += 1
                filtered_vertex_set.append(v)
        else:
            while v.stem_loop_score < 3:
                v.stem_length -= 1
                v.set_stem_loop_score()
            if v.distance >= dist_min and v.distance <= dist_max and v.stem_loop_score >= sl_min and v.stem_loop_score <= sl_max:
                v.reset_info()
                v.set_index(ct)
                ct += 1
                filtered_vertex_set.append(v)

    return filtered_vertex_set


def get_top_alignments(l: int, top_alignment_inds: list, vertex_set: list,
                       cliques: list) -> list:
    """
    Get the alignments from the clique with maximum energy.

    Args:
        l (int): length of the given RNA sequence.
        top_alignment_inds (list): a list of index of cliques that has maximum energy.
        vertex_set (list): a set of vertex.
        cliques (list): a list of possible cliques.

    Returns:
        list: a list of string alignments.
    """
    alignments = []
    # to distinguish different stems, different signs are used to represent predicted alignments
    pair_signs = [('(', ')') , ('[', ']'), ('{', '}')]
    for i in top_alignment_inds:
        # take one top alignment (there might be multiple sinces some alignments may have the same energy)
        alignment = ['.'] * l
        c = cliques[i]
        for j in range(len(c)):
            v = vertex_set[c[j]]
            for k in range(v.stem_length):
                alignment[v.info[k][0]] = pair_signs[j % 3][0]
                alignment[v.info[k][1]] = pair_signs[j % 3][1]

        alignments.append(''.join(alignment))

    return alignments


def bron_kerbosch(a: np.ndarray) -> np.ndarray:
    """
    Find maximal cliques using Bron-Kerbosch algorithm
    The Bron Kerbosch algorithm is an algorithm for finding maximal cliques in an undirected graph.
    
    Reference: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
    of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
    pp: 575?577, September 1973.

    Reference: Cazals, F. and Karande, C., "A note on the problem of reporting 
    maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
    no. 1-3, pp: 564-568, November 2008.
    
    Args:
        a (np.ndarray): an adjacency matrix.

    Returns:
        np.ndarray: a matrix that contains the maximal cliques in its 
    columns.
    """

    n = a.shape[0]
    MC = []  # a list of maximal cliques
    R = np.array([])  # store the currently growing cliques
    P = np.arange(n)  # pro spective nodes connected to all nodes in R
    X = np.array([])  # nodes already processed

    def BK(R, P, X):
        if P.shape[0] == 0 and X.shape[0] == 0:
            newMC = np.zeros(n)
            newMC[R] = 1
            MC.append(newMC)
        else:
            ppivots = np.union1d(P, X).astype(int)
            binP = np.zeros(n)
            binP[P] = 1
            pcounts = a[ppivots, :] @ binP
            ind = np.argmax(pcounts)
            v_p = ppivots[ind]
            temp = np.intersect1d(np.where(a[v_p, :] == 0), P)
            for j in range(len(temp)):
                v = temp[j]
                P = np.setxor1d(P, v)
                Rnew = np.append(R, v)
                Nv = np.where(a[v, :])
                Pnew = np.intersect1d(P, Nv)
                Xnew = np.intersect1d(X, Nv)
                BK(Rnew.astype(int), Pnew, Xnew)
                X = np.block([X, v])
        return

    BK(R, P, X)
    MC = np.array(MC).transpose((1, 0))
    return MC


def compute_energy(cliques, v_set):
    energy_list = []
    for clique in cliques:
        sum_length = 0
        for i in clique:
            sum_length += v_set[i].stem_length
        energy_list.append(sum_length)
    energy_list = np.array(energy_list)
    answer = np.where(energy_list == np.max(energy_list))[0]
    return energy_list, answer


# def set_up_paramters(sequence_type):
#     para = {}
#     if sequence_type == 'short RNA':
#         # minimum stem length
#         para['l'] = 3
#         # lower bound of stem-loop score
#         para['sl1'] = 2
#         # upper bound of stem-loop score
#         para['sl2'] = 20
#         # if a sequence allow pseudoknots exist
#         para['p'] = False
#         # if a sequence allow Wobble pairs exist (G-U)
#         para['w'] = False
#         # if a sequence allow U-U pairs exsit
#         para['uu'] = False  # you can customize this

#     elif sequence_type == 'unknown':
#         para['l'] = 3
#         para['sl1'] = -float('inf')
#         para['sl2'] = float('inf')

#     elif sequence_type == 'tRNA':
#         para['l'] = 3
#         para['sl1'] = 3
#         para['sl2'] = 5.4
#         para['d1'] = 12
#         para['d2'] = 18
#         para['p'] = False
#         para['w'] = True

#     elif sequence_type == '5S-rRNA-Archaeal':
#         # other parameters are customized inside vertex finding function
#         para['w'] = True
#         para['p'] = False
#         para['vertex structure'] = set_parameters_5s_rRNA_archaeal()

#     elif sequence_type == '5S-rRNA-Bacterial':
#         # other parameters are customized inside vertex finding function
#         para['w'] = True
#         para['p'] = False
#         para['vertex structure'] = set_parameters_5s_rRNA_bacterial()

#     return para
