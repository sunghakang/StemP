import numpy as np
from utils.helpers import is_base_pair, bron_kerbosch
from utils.helpers import filter_based_on_stem_loop_score
from utils.helpers import filter_based_on_distance_and_open_ended_modification

class Stem:

    def __init__(self, info=[]):
        """
        Initialize a Stem object.
        """
        self.info = info
        self.neighbors = []
        self.id = 0

    def set_info(self, v: list = []):
        """
        Set up the detailed index of paired bases in a given stem.

        Args:
            v (list, optional): a list of paired bases. For example, [(1,10), 
            (2,9), (3,8)]. Defaults to [].
        """
        self.info = v
        return

    def reset_info(self):
        """
        Reset the detailed index of paired bases in a given stem in terms of the l
        ength of this stem. This function will be used during partial stem modification.
        """
        self.info = self.info[:self.stem_length]
        return

    def set_stem_length(self):
        """
        Set up the stem-length of the stem.
        """
        self.stem_length = len(self.info)
        return

    def set_distance(self, d):
        """
        Set up the distance of a given stem manually.
        """
        self.distance = d
        return

    def set_stem_loop_score(self):
        """
        Calculate and set up the stem-loop-score of a given stem.
        """
        self.stem_loop_score = self.distance / self.stem_length
        return

    def set_index(self, i):
        self.id = i
        return

    def add_neighbor(self, neighbor: list):
        """
        Add neighboring vertex/stems to an existing stem.

        Args:
            neighbor (list): a list of Stem objects.
        """
        if neighbor.id not in self.neighbors:
            self.neighbors.append(neighbor.id)
            neighbor.neighbors.append(self.id)

    def print_info(self):
        """
        Print detailed information about a Stem.
        """
        print('Info: ', self.info, ', d: ', self.distance, ', L: ',
              self.stem_length, ', sl: ', self.stem_loop_score, ', id: ',
              self.id)
        return


class Spectrum_Graph:

    def __init__(self):
        """
        Initialize a Spectrum Graph.
        """
        self.vertices = {}
        self.num_vertices = 0
        self.edges = []
        self.possible_vertex_set = []

    def __iter__(self):
        return iter(self.vertices.values())

    def add_vertex(self, stem: Stem):
        """
        Add a vertex to this graph. Note: each vertex should have a unique id.

        Args:
            stem (Stem): a Stem object.
        """
        self.vertices[stem.id] = []
        self.num_vertices += 1
        return

    def add_vertices(self):
        """
        Add all the vertices/stems from possible_vertex_set into the Spectrum Graph.
        """
        for s in self.possible_vertex_set:
            self.add_vertex(s)

    def add_edge(self, frm: Stem, to: Stem):
        """
        Add an undirected edge between two Stems. 

        Args:
            frm (Stem): The first stem.
            to (Stem): The second stem.
        """
        frm.add_neighbor(to)
        self.vertices[frm.id] = frm.neighbors
        self.vertices[to.id] = to.neighbors
        return

    def add_edges(self, edges: list):
        """
        Add multiple edges from a list into the Spectrum Graph.

        Args:
            edges (list): A list of paired Stems.
        """
        for edge in edges:
            self.add_edge(edge[0], edge[1])

    def set_adjacencyList(self) -> dict:
        """
        Set up the adjacency dictionary list for vertices in the Spectrum Graph.

        Returns:
            dict: a list of adjacency vertices in the Spectrum Graph.
        """
        if len(self.vertices) >= 1:
            return [
                str(key) + ":" + str(self.vertices[key])
                for key in self.vertices.keys()
            ]
        else:
            return dict()

    def set_adjacencyMatrix(self) -> np.ndarray:
        """
        Set up an adjacency matrix for vertices in the Spectrum Graph.

        Returns:
            np.ndarray: an adjacency matrix.
        """
        if len(self.vertices) >= 1:
            self.vertex_names = sorted(self.vertices.keys())
            self.vertex_indices = dict(
                zip(self.vertex_names, range(len(self.vertex_names))))

            self.adjacency_matrix = np.zeros(shape=(len(self.vertices),
                                                    len(self.vertices)))

            for i in range(len(self.vertex_names)):
                for j in range(i, len(self.vertices)):
                    for el in self.vertices[self.vertex_names[i]]:
                        j = self.vertex_indices[el]
                        self.adjacency_matrix[i, j] = 1
            return self.adjacency_matrix
        else:
            return np.array([])

    def save_cliques(self, MC: np.ndarray):
        """
        Find cliques from using MC and save them into self.cliques.
        """
        self.cliques = []
        row, column = MC.shape
        for i in range(column):
            x = []
            for j in range(row):
                if MC[j, i] == 1:
                    x = x + [j]
            self.cliques.append(x)
        return

    def build_cliques(self):
        """
        Find and build cliques from the vertices in the spectrum graph.
        """
        if len(self.edges) == 0:
            self.cliques = [[i] for i in range(len(self.possible_vertex_set))]
        else:
            self.add_edges(self.edges)
            a = self.set_adjacencyMatrix()
            MC = bron_kerbosch(a)
            self.save_cliques(MC)
        return

    def find_stems(self, rna: str, seq_type: str, para: dict):
        """
        Find stems/vertices of a given RNA sequence.

        Args:
            rna (str): an RNA sequence.
            seq_type (str): sequence type.
            para (dict): parameters used in prediction.
        """
        if seq_type == 'short RNA':
            self.find_vertex_v1(rna, para['l'], para['w'], para['uu'])
            self.possible_vertex_set = filter_based_on_stem_loop_score(
                self.possible_vertex_set, para['sl1'], para['sl2'])

        elif seq_type == 'unknown':
            self.find_vertex_v1(rna, para['l'], False, False)
            self.possible_vertex_set = filter_based_on_stem_loop_score(
                self.possible_vertex_set, para['sl1'], para['sl2'])

        elif seq_type == 'tRNA':
            self.find_vertex_v1(rna, para['l'], para['w'], False)
            self.possible_vertex_set = filter_based_on_distance_and_open_ended_modification(
                self.possible_vertex_set, para['sl1'], para['sl2'], para['d1'],
                para['d2'], len(rna))

        elif seq_type in ['5S-rRNA-Archaeal', '5S-rNRA-Bacterial']:
            self.possible_vertex_set = []
            vs = para['vertex_structure']
            b1, b2_1, b2_2, b3_1, b3_2, b2, b3 = vs['b1'], vs['b2_1'], vs[
                'b2_2'], vs['b3_1'], vs['b3_2'], vs['b2'], vs['b3']
            b1_vertex = self.find_structured_vertex(rna, para, b1, 'b1')
            b2_vertex = self.construct_branches(rna, para, b2_1, b2_2, b2,
                                                'b2')
            b3_vertex = self.construct_branches(rna, para, b3_1, b3_2, b3,
                                                'b3')
            self.possible_vertex_set = b1_vertex + b2_vertex + b3_vertex
            self.remove_duplicate_vertex()
            self.reset_vertex_index()

        return

    def construct_branches(self, rna: str, para: dict, para_set_1: dict,
                           para_set_2: dict, parameter: list,
                           label: str) -> list:
        """
        Construct two sets of branches based on given parameters and save all 
        possible branches into a list. Each branch is now a single vertex.

        Args:
            rna (str): an RNA sequence.
            para (dict): parameters used in prediction.
            para_set_1 (dict): parameters of vertex in set 1.
            para_set_2 (dict): parameters of vertex in set 2.
            parameter (list): parameter of a branch.
            label (str): label of a branch.

        Returns:
            list: a list of branch(new vertex) after merging.
        """
        set_1 = self.find_structured_vertex(rna, para, para_set_1)
        set_2 = self.find_structured_vertex(rna, para, para_set_2)
        set_of_branch = self.connect_two_vertex(set_1, set_2, parameter, label)
        return set_of_branch

    def remove_duplicate_vertex(self):
        """
        Remove the duplicated vertex in possible_vertex_set.
        """
        vertex_info_set = set()
        repeat_ind = []
        for i in range(len(self.possible_vertex_set)):
            v = self.possible_vertex_set[i]
            if str(v.info) not in vertex_info_set:
                vertex_info_set.add(str(v.info) + str(v.label))
            else:
                repeat_ind = [i] + repeat_ind
        for i in repeat_ind:
            del self.possible_vertex_set[i]
        return

    def reset_vertex_index(self):
        """
        Reset the id for each vertex in possible_vertex_set.
        """
        count = 0
        for v in self.possible_vertex_set:
            v.id = count
            count += 1
        return

    def connect_two_vertex(self, v1_set: list, v2_set: list, branch_para: dict,
                           label: str) -> list:
        """
        Connect one vertice from v1_set and another vertice from v2_set to form a branch. 
        All the branches are saved in vertex_set.

        Args:
            v1_set (list): a list of vertex/stems. (The first part of a branch.)
            v2_set (list): another list of vertex/stem. (The second part of a branch.)
            branch_para (dict): parameters of a branch to be formed.
            label (str): a uniform label for the branch.

        Returns:
            list: a list of branches. Each element in this list is considered as a new vertex.
        """
        sl_min, sl_max = branch_para['sl_min'], branch_para['sl_max']
        vertex_set = []
        ct = [0]
        for i in range(len(v1_set)):
            v1 = v1_set[i]
            for j in range(len(v2_set)):
                v2 = v2_set[j]
                if v1.info[-1][0] < v2.info[0][0] and v1.info[-1][1] > v2.info[
                        0][1]:
                    new_vertex_info = v1.info + v2.info
                    self.check_and_save_vertex(vertex_set, sl_min, sl_max,
                                               -float('inf'), float('inf'), ct,
                                               new_vertex_info.copy(), label)
        return vertex_set

    def find_structured_vertex(self,
                               rna: str,
                               para: dict,
                               branch_para: dict,
                               label=None) -> list:
        """
        Find some vertex that satisfies a specific structure specified in branch_para.

        Args:
            rna (str): an RNA sequence.
            para (dict): general parameter for a vertex.
            branch_para (dict): parameter for a branch.
            label (_type_, optional): label for a branch. Defaults to None.

        Returns:
            list: a list of branches
        """
        branch_vertex = []
        for key, vertex_structures in branch_para.items():
            if key == 'type-1 vertex':
                for vertex_structure in vertex_structures:
                    branch_vertex += self.find_one_piece_vertex(
                        rna, para['w'], vertex_structure, label)
            elif key == 'type-2 vertex':
                for vertex_structure in vertex_structures:
                    branch_vertex += self.find_two_piece_vertex(
                        rna, para['w'], vertex_structure, label)
            elif key == 'type-3 vertex':
                for vertex_structure in vertex_structures:
                    branch_vertex += self.find_three_piece_vertex(
                        rna, para['w'], vertex_structure, label)
        return branch_vertex

    def find_vertex_v1(self, rna: str, stack_length: int, has_wobble: bool,
                       has_uu: bool):
        """
        Find some vertex that satisfies certain conditions.

        Args:
            rna (str): an RNA sequence.
            stack_length (int): minimum length of stem/vertex.
            has_wobble (bool): whether wobble pair G-U is allowed.
            has_uu (bool): whether U-U pair is allowed.
        """
        for i in range(len(rna)):
            for j in range(i + stack_length, len(rna)):
                starting, ending, new_vertex_info = i, j, []
                while is_base_pair(rna[starting], rna[ending], has_wobble,
                                   has_uu) and starting < ending:
                    x = (starting, ending)
                    new_vertex_info.append(x)
                    starting += 1
                    ending -= 1

                if len(new_vertex_info) >= stack_length:
                    new_vertex = Stem(new_vertex_info)
                    new_vertex.set_distance(new_vertex_info[0][1] -
                                            new_vertex_info[0][0])
                    new_vertex.set_stem_length()
                    new_vertex.set_stem_loop_score()
                    self.possible_vertex_set.append(new_vertex)

        return

    def find_one_piece_vertex(self, rna: str, has_wobble: bool,
                              vertex_structure: dict, label) -> list:
        """
        Find a vertex with consecutive base pairs.

        Args:
            rna (str): an RNA sequence.
            has_wobble (bool): whether wobble pair G-U is allowed.
            vertex_structure (dict): a dictionary of parameters about the structure of verteices.
            label (None or str): label of a vertex.

        Returns:
            list: a list of vertex.
        """
        l_min, l_max, sl_min, sl_max = vertex_structure[
            'l_min'], vertex_structure['l_max'], vertex_structure[
                'sl_min'], vertex_structure['sl_max']
        vertex_set = []
        ct = [0]
        for i in range(len(rna)):
            for j in range(i + l_min, len(rna)):
                starting, ending, new_vertex_info = i, j, []
                while is_base_pair(
                        rna[starting], rna[ending], has_wobble, False
                ) and starting < ending and len(new_vertex_info) < l_max:
                    x = (starting, ending)
                    new_vertex_info.append(x)
                    starting += 1
                    ending -= 1
                    self.check_and_save_vertex(vertex_set, sl_min, sl_max,
                                               l_min, l_max, ct,
                                               new_vertex_info.copy(), label)
        return vertex_set

    def find_two_piece_vertex(self, rna: str, has_wobble: bool,
                              vertex_structure: dict, label) -> list:
        """
        Find a vertex with two disconnected sub-stems.

        Args:
            rna (str): an RNA sequence.
            has_wobble (bool): whether a wobble pair G-U is allowed.
            vertex_structure (dict): a dictionary of parameters about the structure of vertices.
            label (None or str): label of a vertex.

        Returns:
            list: a list of vertex.
        """
        l1, l2, sl_min, sl_max, gap_1, gap_2 = vertex_structure[
            'l1'], vertex_structure['l2'], vertex_structure[
                'sl_min'], vertex_structure['sl_max'], vertex_structure[
                    'gap_1'], vertex_structure['gap_2']
        vertex_set = []
        total_length = l1 + l2
        stack_length = 2
        ct = [0]
        for i in range(len(rna)):
            for j in range(i + stack_length, len(rna)):
                starting, ending, new_vertex_info = i, j, []
                k = 0
                while is_base_pair(rna[starting], rna[ending], has_wobble,
                                   False) and starting < ending and len(
                                       new_vertex_info) < total_length:
                    k += 1
                    x = (starting, ending)
                    new_vertex_info.append(x)
                    if k != l1:
                        starting += 1
                        ending -= 1
                    else:
                        for k1 in range(gap_1):
                            for k2 in range(gap_2):
                                if ending - (k2 + 1) >= 0 and (is_base_pair(
                                        rna[starting + k1 + 1], rna[ending -
                                                                    (k2 + 1)],
                                        has_wobble, False) or is_base_pair(
                                            rna[starting + k1 + 1],
                                            rna[ending -
                                                (k2 + 1)], has_wobble, False)):
                                    continue
                        starting += 1 + gap_1
                        ending -= (1 + gap_2)
                    if starting > len(rna) - 1 or ending <= starting:
                        break
                self.check_and_save_vertex(vertex_set, sl_min, sl_max,
                                           total_length, total_length, ct,
                                           new_vertex_info.copy(), label)

        return vertex_set

    def find_three_piece_vertex(self, rna: str, has_wobble: bool,
                                vertex_structure: dict, label) -> list:
        """
        Find a vertex with three disconnected sub-stems.

        Args:
            rna (str): an RNA sequence.
            has_wobble (bool): whether wobble pair G-U is allowed.
            vertex_structure (dict): a dictionary of parameters about the structure of verteices.
            label (None or str): label of a vertex.

        Returns:
            list: a list of vertex.
        """
        l1, l2, l3, sl_min, sl_max, gap_id_1, gap_id_2, gap_1, gap_2 = vertex_structure[
            'l1'], vertex_structure['l2'], vertex_structure[
                'l3'], vertex_structure['sl_min'], vertex_structure[
                    'sl_max'], vertex_structure['gap_id_1'], vertex_structure[
                        'gap_id_2'], vertex_structure[
                            'gap_1'], vertex_structure['gap_2']
        total_length = l1 + l2 + l3
        vertex_set = []
        ct = [0]
        for i in range(len(rna)):
            for j in range(i + 2, len(rna)):
                starting, ending, new_vertex_info = i, j, []
                k = 0
                while (ending >= 0 and is_base_pair(
                        rna[starting], rna[ending], has_wobble,
                        False)) and (starting < ending
                                     and len(new_vertex_info) < total_length):
                    k += 1
                    x = (starting, ending)
                    new_vertex_info.append(x)
                    if k == gap_id_1:
                        starting += (1 + gap_1[0])
                        ending -= (1 + gap_1[1])
                    elif k == gap_id_2:
                        starting += (1 + gap_2[0])
                        ending -= (1 + gap_2[1])
                    else:
                        starting += 1
                        ending -= 1
                self.check_and_save_vertex(vertex_set, sl_min, sl_max,
                                           total_length, total_length, ct,
                                           new_vertex_info.copy(), label)
        return vertex_set

    def check_and_save_vertex(self,
                              vertex_set: list,
                              sl_min: float,
                              sl_max: float,
                              l_min: int,
                              l_max: int,
                              ct: list,
                              new_vertex_info: list,
                              label=None):
        """
        Check if a vertex satisfies a specific condition and save it to  vertex_set.

        Args:
            vertex_set (list): a list of candidate vertex.
            sl_min (float): minimum stem-loop-score.
            sl_max (float): maximum stem-loop-score.
            l_min (int): minimum stem length.
            l_max (int): maximum stem length.
            ct (list): current total count of possible vertex.
            new_vertex_info (list): detailed vertex base information of the current vertex.
            label (_type_, optional): label of a vertex. Defaults to None.
        """
        if len(new_vertex_info) >= l_min and len(new_vertex_info) <= l_max:
            new_vertex = Stem(new_vertex_info)
            new_vertex.set_distance(new_vertex_info[0][1] -
                                    new_vertex_info[0][0])
            new_vertex.set_stem_length()
            new_vertex.set_stem_loop_score()
            if new_vertex.stem_loop_score >= sl_min and new_vertex.stem_loop_score <= sl_max:
                new_vertex.id = ct[0]
                ct[0] += 1
                if label is not None:
                    new_vertex.label = label
                vertex_set.append(new_vertex)
        return

    def find_edge(self, has_pseu: bool, has_three_branches=False):
        """
        Find edges between the vertices in the spectrum graph.

        Args:
            has_pseu (bool): whether pseudknots are allowed.
            has_three_branches (bool, optional): whether the structure has only 3 branches. Defaults to False.
        """
        for i in range(len(self.possible_vertex_set)):
            i_info = self.possible_vertex_set[i].info
            for j in range(i + 1, len(self.possible_vertex_set)):
                if not (has_three_branches is True
                        and self.possible_vertex_set[i].label
                        == self.possible_vertex_set[j].label):
                    j_info = self.possible_vertex_set[j].info
                    if has_pseu is False:
                        if i_info[0][1] < j_info[0][0] or j_info[0][
                                1] < i_info[0][0] or (
                                    i_info[-1][0] < j_info[0][0]
                                    and j_info[0][1] < i_info[-1][1]) or (
                                        j_info[-1][0] < i_info[0][0]
                                        and i_info[0][1] < j_info[-1][1]):
                            s = (self.possible_vertex_set[i],
                                 self.possible_vertex_set[j])
                            self.edges.append(s)
                    else:
                        if i_info[0][1] < j_info[0][0] or (
                                i_info[-1][0] < j_info[0][0]
                                and j_info[0][1] < i_info[-1][1]) or (
                                    i_info[-1][0] < j_info[0][0]
                                    and j_info[-1][0] < i_info[-1][1]
                                    and i_info[0][1] < j_info[-1][1]):
                            s = (self.possible_vertex_set[i],
                                 self.possible_vertex_set[j])
                            self.edges.append(s)
        return