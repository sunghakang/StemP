from utils.helpers import get_top_alignments, compute_energy
from spectrum_graph import Spectrum_Graph
import time


def stemP_rna_structure_pred(rna: str, seq_type: str, para: dict) -> list:
    """This function performs prediction of RNA secondary structure given an RNA sequence.

    Args:
        rna (str): an RNA sequence.
        seq_type (str): sequence type.
        para (dict): parameters used in prediction.

    Returns:
        list: a list of possible alignments with top energy.
    """
    start_time = time.time()

    # build spectrum graph
    graph = Spectrum_Graph()
    graph.find_stems(rna, seq_type, para)
    graph.add_vertices()
    print('Number of vertices found : ', len(graph.possible_vertex_set))
    graph.find_edge(para['p'], seq_type == '5S-rRNA-Archaeal')
    print('Number of edges found  : ', len(graph.edges))

    graph.build_cliques()
    print('Number of cliques found: ', len(graph.cliques))

    # compute energy for each clique
    _, top_alignment_inds = compute_energy(graph.cliques,
                                           graph.possible_vertex_set)

    # find the alignments with top energy
    top_alignments = get_top_alignments(len(rna), top_alignment_inds,
                                        graph.possible_vertex_set,
                                        graph.cliques)

    run_time = time.time() - start_time
    print('Running time: ', run_time, ' seconds.')

    return top_alignments