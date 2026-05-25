## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Contact: lipshitz@math.columbia.edu, sucharit@math.columbia.edu
## Pure Python port (C) 2026.

## This file is part of KhovanovSteenrod.

## KhovanovSteenrod is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## KhovanovSteenrod is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with KhovanovSteenrod; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

import numpy as np
from copy import deepcopy
import time

from gf2 import GF2Matrix, _pivot_rows, _left_null_space, _row_space, _rank


################## Some easy (but convenient) functions. ######################

def first_difference(a, b):
    "Returns index of first place tuples a and b differ."
    if a[0] != b[0]:
        return 0
    return first_difference(a[1:], b[1:]) + 1

def second_difference(a, b):
    "Returns index of second place tuples a and b differ."
    first = first_difference(a, b)
    return first_difference(a[first+1:], b[first+1:]) + first + 1

def second_index(data, elt):
    "Returns the second place that elt occurs in data"
    return data[data.index(elt)+1:].index(elt) + data.index(elt) + 1

def doubletest_keys(a, b, dic, p=0):
    """dic is a dictionary. If p=0, checks if a is in dic, and b is in dic[a]
    (if p=1, checks if a in dic and b in dic[a][1])"""
    if a in dic.keys():
        if p == 0 and b in dic[a]:
            return True
        if p == 1 and b in dic[a][1]:
            return True
    return False

def occurs_twice(data, elt):
    # Checks if +-elt occurs twice in a list.
    if len([x for x in data if (x == elt or x == -elt)]) > 1:
        return True
    return False

def hypercube(n, faces=False):
    """
    Inputs: n - integer - dimension of the hypercube
    Returns: (V, E) where V is a list of vertices and E is a list of edges,
    or (V, E, F) if faces=True.
    """
    if n == 1 and faces:
        return ([(0,), (1,)], [((0,), (1,))], [])
    if n == 1 and not faces:
        return ([(0,), (1,)], [((0,), (1,))])
    if n == 2 and faces:
        return (
            [(0, 0), (1, 0), (0, 1), (1, 1)],
            [((0, 0), (1, 0)), ((0, 1), (1, 1)), ((0, 0), (0, 1)), ((1, 0), (1, 1))],
            [((0, 0), (1, 0), (0, 1), (1, 1)),]
        )
    if faces:
        (nMinusOneVerts, nMinusOneEdges, nMinusOneFaces) = hypercube(n-1, True)
    else:
        (nMinusOneVerts, nMinusOneEdges) = hypercube(n-1, False)
    vertices = [vert+(0,) for vert in nMinusOneVerts] + [vert+(1,) for vert in nMinusOneVerts]
    edges = ([(edge[0]+(0,), edge[1]+(0,)) for edge in nMinusOneEdges] +
             [(edge[0]+(1,), edge[1]+(1,)) for edge in nMinusOneEdges] +
             [(vertex+(0,), vertex+(1,)) for vertex in nMinusOneVerts])
    if not faces:
        return (vertices, edges)
    else:
        faces = (
            [(face[0]+(0,), face[1]+(0,), face[2]+(0,), face[3]+(0,)) for face in nMinusOneFaces] +
            [(face[0]+(1,), face[1]+(1,), face[2]+(1,), face[3]+(1,)) for face in nMinusOneFaces] +
            [(edge[0]+(0,), edge[1]+(0,), edge[0]+(1,), edge[1]+(1,)) for edge in nMinusOneEdges]
        )
        return (vertices, edges, faces)

def vertices_out(v):
    "Returns a list of the vertices at the ends of outgoing edges from the vertex v."
    return [v[0:i]+(1,)+v[i+1:] for i in range(len(v)) if v[i] == 0]

def vertices_2out(v):
    "Returns a list of triples (m1,m2,f) for all outgoing 2D faces from v."
    ans = []
    for i in range(len(v)):
        if v[i] == 0:
            m1 = v[0:i]+(1,)+v[i+1:]
            for j in range(i+1, len(v)):
                if v[j] == 0:
                    m2 = v[0:j]+(1,)+v[j+1:]
                    f  = m1[0:j]+(1,)+m1[j+1:]
                    ans.append((m1, m2, f))
    return ans

###############################################################################


############## The class CycList and some related functions. ##################

class CycList():
    """
    Input:
    data: a list or tuple of integers.
    """

    def __init__(self, data):
        self.data = list(data)

    def __hash__(self):
        return hash(tuple(self.data))

    def __eq__(self, other):
        return self.data == other.data

    def __repr__(self):
        return "CycList(" + repr(self.data) + ")"

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[item]

    def rotate(self, n):
        "Return result of rotating leftwards by n units."
        return CycList(self.data[n:] + self.data[0:n])

    def reflect(self):
        n = len(self.data)
        return CycList([-self.data[n-i-1] for i in range(len(self.data))])

    def equiv(self, other):
        md = list(self.data)
        od = list(other.data)
        nmd = list(self.reflect().data)
        md.sort()
        od.sort()
        nmd.sort()
        if md != od and nmd != od:
            return False
        else:
            md = list(self.data)
            if nmd == od:
                md = list(self.reflect().data)
            od = list(other.data)
            if od[(od.index(md[0])+1) % (len(od))] == md[1 % (len(od))]:
                return True
            else:
                return False

    def merge(self, other, n):
        "Merges cyclic lists self and other along the edge labeled n. Returns: CycList"
        if not (n in self.data):
            new_a = self.reflect()
        else:
            new_a = CycList(self.data)
        if not (n in other.data):
            new_b = other.reflect()
        else:
            new_b = CycList(other.data)
        a_loc = list(new_a.data).index(n)
        b_loc = list(new_b.data).index(n)
        answer = (new_a.data[a_loc+1:] + new_a.data[0:a_loc] + [-n,] +
                  new_b.data[b_loc+1:] + new_b.data[0:b_loc] + [-n,])
        return CycList(answer)

    def split(self, n):
        "Splits cyclic list self along the edge labeled n. Returns (b,c) of CycList's"
        adata = list(self.data)
        newn = n
        if not (n in adata):
            newn = -n
        first  = adata.index(newn)
        second = first + adata[first+1:].index(newn) + 1
        b = CycList(adata[first+1:second] + [-newn,])
        c = CycList(adata[second+1:] + adata[0:first] + [-newn,])
        return (b, c)


def equiv(a, b):
    return a.equiv(b)

def merge(a, b, n):
    return a.merge(b, n)

def split(a, n):
    return a.split(n)


def ladybug(start, mid, edges):
    """
    start is a CycList, mid is a list of 4 CycLists, edges is a list of two
    edge labels whose endpoints are linked S^0s on start (ladybug configuration).
    Returns (matching_dict, correction) where matching_dict is a self-matching
    of {0,1,2,3} and correction is 0 or 1.
    """
    answer = dict()
    a = abs(edges[0])
    b = abs(edges[1])

    # Rotate/reflect so that start looks like (a, ?, -b, ?, a, ?, -b, ?)
    if a in start.data:
        tempns = CycList(start.data)
    else:
        tempns = start.reflect()
    ns = tempns.rotate(tempns.data.index(a))

    # Adjust mids so that a (positive) occurs in all mids as first entry
    nmids = []
    for x in mid:
        if a in x.data:
            nmids.append(x.rotate(x.data.index(a)))
        else:
            temp = x.reflect()
            nmids.append(temp.rotate(temp.data.index(a)))

    # Case: one of the positively oriented intervals is nonempty
    if ns.data[1] != -b or ns.data[second_index(ns.data, a)+1] != -b:
        if ns.data[1] == -b:
            ns = ns.rotate(second_index(ns.data, a))
        c = ns.data[1]
        matched_pair = [x for x in nmids if (x.data[1] == c or x.data[len(x.data)-1] == -c)]
        answer[nmids.index(matched_pair[0])] = nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])] = nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if i not in answer.keys()]
        answer[unmatched[0]] = unmatched[1]
        answer[unmatched[1]] = unmatched[0]
        return (answer, 0)

    # Case: both positively oriented intervals are empty
    else:
        if second_index(ns.data, a) == 2:
            ns = ns.rotate(second_index(ns.data, a))
        c = ns.data[2]
        matched_pair = [
            x for x in nmids
            if (x.data[(x.data.index(b)+1) % len(x.data)] == c or
                x.data[len(x.data)-2] == -c)
        ]
        answer[nmids.index(matched_pair[0])] = nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])] = nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if i not in answer.keys()]
        answer[unmatched[0]] = unmatched[1]
        answer[unmatched[1]] = unmatched[0]
        return (answer, 1)

###############################################################################


############### The class ResolvedKnot and related function. ##################

class ResolvedKnot():
    """
    Inputs:
     --vertex_incidences: a list of CycList's

    Data:
    N = number of edges (= number of crossings in unresolved knot)
    C = number of components
    """

    def __init__(self, vertex_incidences, label=[]):
        self.vertex_incidences = list(vertex_incidences)
        self.N = len(sum([list(x.data) for x in vertex_incidences], [])) // 2
        if label:
            self.label = tuple(label)
        else:
            self.label = self.N * (0,)
        self.edge_incidences = dict([(i+1, list()) for i in range(self.N)])
        for x in vertex_incidences:
            for y in x.data:
                self.edge_incidences[abs(y)].append(vertex_incidences.index(x))
        self.C = len(self.vertex_incidences)

    def __repr__(self):
        return ("ResolvedKnot with label: " + repr(self.label) +
                "\n     v.i.: " + repr(self.vertex_incidences) + "\n")

    def flip(self, edge):
        "Returns the resolved knot obtained by switching the resolution at edge."
        new_label = self.label[0:edge-1] + ((self.label[edge-1]+1) % 2,) + self.label[edge:]
        if self.edge_incidences[edge][0] == self.edge_incidences[edge][1]:
            # Split
            sc = self.edge_incidences[edge][0]
            other_components = self.vertex_incidences[0:sc] + self.vertex_incidences[sc+1:]
            new_vertex_incidences = tuple(other_components) + split(self.vertex_incidences[sc], edge)
        else:
            # Merge
            mc1 = min(self.edge_incidences[edge][0], self.edge_incidences[edge][1])
            mc2 = max(self.edge_incidences[edge][0], self.edge_incidences[edge][1])
            other_components = (self.vertex_incidences[0:mc1] +
                                self.vertex_incidences[mc1+1:mc2] +
                                self.vertex_incidences[mc2+1:])
            new_vertex_incidences = (tuple(other_components) +
                                     (merge(self.vertex_incidences[mc1],
                                            self.vertex_incidences[mc2], edge),))
        return ResolvedKnot(new_vertex_incidences, new_label)

###############################################################################


############ Now some main functions. #########################################

def populate_vertices(zero_res):
    "Returns: dictionary with keys elements of hypercube(N)[0] and values ResolvedKnot"
    N = zero_res.N
    cube = dict()
    cube[N*(0,)] = zero_res
    for s in hypercube(N)[0][1:]:
        i = list(s).index(1)
        prev_s = s[0:i]+(0,)+s[i+1:]
        cube[s] = cube[prev_s].flip(i+1)
    return cube

def find_fly_in_soup(fly, soup):
    "Returns the index of the CycList fly in the list of CycList's soup"
    for b in soup:
        if equiv(fly, b):
            return list(soup).index(b)
    raise Exception("Can't find fly " + repr(fly) + " in soup " + repr(soup))

def populate_edges(cube):
    """
    Input: output of populate_vertices().
    Returns: dictionary keyed by edges of N-D hypercube with values
    ('m'/'s', a, b, c, d, bijection).
    """
    N = len(next(iter(cube.keys())))
    answer = dict()
    for e in hypercube(N)[1]:
        initial_vert = cube[e[0]]
        final_vert   = cube[e[1]]
        edge_flipping = first_difference(e[0], e[1]) + 1
        if initial_vert.C > final_vert.C:
            move_type = 'm'
        else:
            move_type = 's'
        a = min(initial_vert.edge_incidences[edge_flipping][0],
                initial_vert.edge_incidences[edge_flipping][1])
        b = max(initial_vert.edge_incidences[edge_flipping][0],
                initial_vert.edge_incidences[edge_flipping][1])
        c = min(final_vert.edge_incidences[edge_flipping][0],
                final_vert.edge_incidences[edge_flipping][1])
        d = max(final_vert.edge_incidences[edge_flipping][0],
                final_vert.edge_incidences[edge_flipping][1])
        bijection = dict()
        for x in range(initial_vert.C):
            if x != a and x != b:
                bijection[x] = find_fly_in_soup(initial_vert.vertex_incidences[x],
                                                final_vert.vertex_incidences)
        answer[e] = (move_type, a, b, c, d, bijection)
    return answer

def chain_groups(cube, t_q_correct=(0, 0)):
    """
    Returns dictionary keyed by bigrading (t,q) with values lists of generators.
    """
    (t_correct, q_correct) = t_q_correct
    generators = dict()
    for s in cube.keys():
        for g in hypercube(cube[s].C)[0]:
            t = sum(s) + t_correct
            q = len(g) - 2*sum(g) + t + q_correct - t_correct
            if (t, q) in generators:
                generators[(t, q)].append((s, g))
            else:
                generators[(t, q)] = [(s, g),]
    return generators

def chain_maps(group, edges, print_progress=True):
    """
    Returns dictionary keyed by (t,q) with values GF2Matrix (chain map to (t+1,q)).
    """
    maps = dict()
    group_keys = list(group.keys())
    for idx, k in enumerate(group_keys):
        if print_progress:
            print(repr(idx+1) + '/' + repr(len(group_keys)))
        l = (k[0]+1, k[1])
        if l in group:
            maps[k] = GF2Matrix(len(group[k]), len(group[l]))
        for i in range(len(group[k])):
            (s, x) = group[k][i]
            for t in vertices_out(s):
                (edge_type, a, b, c, d, bijection) = edges[(s, t)]
                if edge_type == 'm':
                    if x[a] != 1 or x[b] != 1:
                        val = (len(x)-1) * [0,]
                        for dummy in list(range(0, a)) + list(range(a+1, b)) + list(range(b+1, len(x))):
                            val[bijection[dummy]] = x[dummy]
                        val[c] = x[a] + x[b]   # 0 or 1 (never 2 since not both 1)
                        j = group[l].index((t, tuple(val)))
                        maps[k][i, j] = 1
                if edge_type == 's':
                    val = (len(x)+1) * [0,]
                    for dummy in list(range(0, a)) + list(range(a+1, len(x))):
                        val[bijection[dummy]] = x[dummy]
                    if x[a] == 1:
                        val[c] = 1
                        val[d] = 1
                        j = group[l].index((t, tuple(val)))
                        maps[k][i, j] = 1
                    if x[a] == 0:
                        val[c] = 0
                        val[d] = 1
                        j = group[l].index((t, tuple(val)))
                        maps[k][i, j] = 1
                        val[c] = 1
                        val[d] = 0
                        j = group[l].index((t, tuple(val)))
                        maps[k][i, j] = 1
    return maps

def complexes(groups, maps, print_progress=True):
    """
    Computes homology representatives.

    Returns (None, homans) where homans[q][t] = (rank, cycle_rep_list).
    cycle_rep_list is a list of numpy uint8 1-D arrays, each a cycle
    representative for a homology generator.
    """
    homans = dict()

    # Find t-range for each q
    temp = dict()
    for (t, q) in groups.keys():
        if q in temp:
            temp[q][0] = min(temp[q][0], t)
            temp[q][1] = max(temp[q][1], t)
        else:
            temp[q] = [t, t]

    temp_keys = list(temp.keys())
    for idx, q in enumerate(temp_keys):
        if print_progress:
            print(repr(idx+1) + "/" + repr(len(temp_keys)))
        start_time = time.time()
        homans[q] = dict()

        for t in range(temp[q][0], temp[q][1]+1):
            if (t, q) not in groups:
                homans[q][t] = (0, [])
                continue

            n = len(groups[(t, q)])

            # Kernel of outgoing differential: {v : v @ maps[(t,q)] = 0}
            if (t, q) in maps:
                ker_vecs = maps[(t, q)].kernel_basis()
            else:
                ker_vecs = [np.eye(n, dtype=np.uint8)[k] for k in range(n)]

            if not ker_vecs:
                homans[q][t] = (0, [])
                continue

            # Image of incoming differential: row space of maps[(t-1,q)]
            if (t-1, q) in maps:
                im_vecs = maps[(t-1, q)].row_space_basis()
            else:
                im_vecs = []

            # Find cycle representatives for ker/im via pivot rows
            if im_vecs:
                all_rows = im_vecs + ker_vecs
                all_arr  = np.array([np.asarray(r, dtype=np.uint8) for r in all_rows],
                                    dtype=np.uint8)
                pivots    = _pivot_rows(all_arr)
                im_rank   = len(im_vecs)
                cycle_reps = [all_rows[p] for p in pivots if p >= im_rank]
            else:
                cycle_reps = list(ker_vecs)

            homans[q][t] = (len(cycle_reps), cycle_reps)

        if print_progress:
            print(repr(q) + "::: time " + repr(time.time() - start_time))

    return (None, homans)

###############################################################################


############### The main KnotKh class #########################################

class KnotKh():
    """
    Main computation class for Khovanov homology and Steenrod squares.

    Input:
        data         -- ResolvedKnot (the 0-resolution)
        nplus_nminus -- tuple (nplus, nminus) of crossing numbers
    """

    def __init__(self, data, nplus_nminus=(0, 0), print_progress=True):
        (nplus, nminus) = nplus_nminus
        self.data = ResolvedKnot(data.vertex_incidences)
        if print_progress:
            print("Populating vertices.")
        self.vertices = populate_vertices(self.data)
        if print_progress:
            print("Populating edges.")
        self.edges = populate_edges(self.vertices)
        if print_progress:
            print("Finding chain groups.")
        self.groups = chain_groups(self.vertices, (-nminus, nplus - 2*nminus))

        start_time = time.time()
        if print_progress:
            print("Finding chain maps.")
        self.maps = chain_maps(self.groups, self.edges, print_progress)
        if print_progress:
            print("Time taken is " + repr(time.time() - start_time))

        start_time = time.time()
        if print_progress:
            print("Generating chain complexes.")
        (_, self.newhom) = complexes(self.groups, self.maps, print_progress)
        if print_progress:
            print("Time taken is " + repr(time.time() - start_time))
            print("Done computing Khovanov homology (over F_2).")

        self.sign1_stored = dict()
        self.sign2_stored = dict()
        self.lb_stored    = dict()

        self.hom         = dict()
        self.hom_support = []

        self.Sq1 = dict()
        self.Sq2 = dict()

    def d2(self, print_output=False):
        "Checks d^2 = 0; returns True iff OK."
        are_we_good = True
        for k in self.maps.keys():
            l = (k[0]+1, k[1])
            if l in self.maps:
                d2 = self.maps[k] @ self.maps[l]
                if d2.rank() != 0:
                    are_we_good = False
                if print_output:
                    print(d2.rank())
        return are_we_good

    def hom_rank(self):
        "Returns list of ((t,q), rank) tuples for all nonzero homology groups."
        ranks = []
        for q in self.newhom.keys():
            for t in self.newhom[q].keys():
                r = self.newhom[q][t][0]
                if r > 0:
                    ranks.append(((t, q), r))
        ranks.sort(key=lambda x: (x[0][1], x[0][0]))
        self.hom = dict(ranks)
        self.hom_support = [(t, q) for ((t, q), r) in ranks]
        return ranks

    def print_hom_rank(self):
        "Prints results of hom_rank() nicely."
        answer = self.hom_rank()
        print_rank(answer)

    def __repr__(self):
        return "KNOT(" + repr(self.data) + ")"

    def sq1(self, tq, n=0, print_progress=True, show_chain=False):
        """
        Computes Sq^1 map on bigrading (t,q).
        If n=0, computes for all basis elements and returns a GF2Matrix.
        If n>0, computes for the n-th basis element and returns a numpy array.
        """
        (t, q) = tq
        if n < 0:
            raise Exception("n must be >= 0")
        if n > self.newhom[q][t][0]:
            raise Exception("Rank is " + repr(self.newhom[q][t][0]) +
                            " so can't find Sq1 of element " + repr(n))
        if n == 0:
            if (t, q) in self.Sq1:
                return self.Sq1[(t, q)]
            r = self.newhom[q][t][0]
            l = []
            for i in range(1, r+1):
                l.append(self.sq1((t, q), i, print_progress=print_progress))
            self.Sq1[(t, q)] = GF2Matrix.from_rows(l) if l else GF2Matrix(0, self.newhom[q].get(t+1, (0,[]))[0])
            return self.Sq1[(t, q)]
        else:
            u = self.newhom[q][t][1][n-1]

            if (t, q) in self.maps and self.newhom[q].get(t+1, (0,[]))[0] > 0:
                v = self.chainsq1(u, (t, q), print_progress=print_progress)

                r_tplus1 = self.newhom[q][t+1][0]
                cycle_reps_tplus1 = self.newhom[q][t+1][1]
                boundary_rows = self.maps[(t, q)].rows()

                all_rows = list(cycle_reps_tplus1) + [np.asarray(r, dtype=np.uint8) for r in boundary_rows]
                temp_mat = GF2Matrix.from_rows(all_rows)
                v_in_hom = temp_mat.solve_left(v)
                v_in_hom = v_in_hom[:r_tplus1]

                if show_chain:
                    return (v, v_in_hom)
                return v_in_hom

            elif show_chain:
                n_chain = len(self.groups[(t+1, q)]) if (t+1, q) in self.groups else 0
                n_hom   = self.newhom[q].get(t+1, (0, []))[0]
                return (np.zeros(n_chain, dtype=np.uint8),
                        np.zeros(n_hom, dtype=np.uint8))
            else:
                n_hom = self.newhom[q].get(t+1, (0, []))[0]
                return np.zeros(n_hom, dtype=np.uint8)

    def sq2(self, tq, n=0, print_progress=True, type='oi', show_chain=False):
        """
        Computes Sq^2 on bigrading (t,q).
        If n=0, computes for all basis elements and returns a GF2Matrix.
        If n>0, computes for the n-th basis element and returns a numpy array.
        """
        (t, q) = tq
        if n < 0:
            raise Exception("n must be >= 0")
        if n > self.newhom[q][t][0]:
            raise Exception("Rank is " + repr(self.newhom[q][t][0]) +
                            " so can't find Sq2 of element " + repr(n))
        if n == 0:
            if (t, q) in self.Sq2:
                return self.Sq2[(t, q)]
            r = self.newhom[q][t][0]
            l = []
            for i in range(1, r+1):
                l.append(self.sq2((t, q), i, print_progress=print_progress, type=type))
            self.Sq2[(t, q)] = GF2Matrix.from_rows(l) if l else GF2Matrix(0, self.newhom[q].get(t+2, (0,[]))[0])
            return self.Sq2[(t, q)]
        else:
            u = self.newhom[q][t][1][n-1]

            t_plus1_in_maps = (t+1, q) in self.maps
            t_plus2_rank    = self.newhom[q].get(t+2, (0, []))[0]

            if (t, q) in self.maps and t_plus1_in_maps and t_plus2_rank > 0:
                v = self.chainsq2(u, (t, q), print_progress=print_progress, type=type)

                r_tplus2 = self.newhom[q][t+2][0]
                cycle_reps_tplus2 = self.newhom[q][t+2][1]
                boundary_rows = self.maps[(t+1, q)].rows()

                all_rows = list(cycle_reps_tplus2) + [np.asarray(r, dtype=np.uint8) for r in boundary_rows]
                temp_mat = GF2Matrix.from_rows(all_rows)
                v_in_hom = temp_mat.solve_left(v)
                v_in_hom = v_in_hom[:r_tplus2]

                if show_chain:
                    return (v, v_in_hom)
                return v_in_hom

            elif show_chain:
                n_chain = len(self.groups[(t+2, q)]) if (t+2, q) in self.groups else 0
                n_hom   = self.newhom[q].get(t+2, (0, []))[0]
                return (np.zeros(n_chain, dtype=np.uint8),
                        np.zeros(n_hom, dtype=np.uint8))
            else:
                n_hom = self.newhom[q].get(t+2, (0, []))[0]
                return np.zeros(n_hom, dtype=np.uint8)

    def chainlady_map(self, u, tq, print_progress=True):
        """
        u is cycle in bigrading (t,q); returns numpy uint8 array (cycle in (t+2,q)).
        """
        (t, q) = tq
        if print_progress:
            print("Finding the ladybugs.")

        if (t, q) in self.lb_stored:
            lb = self.lb_stored[(t, q)]
        else:
            lb = GF2Matrix(len(self.groups[(t, q)]), len(self.groups[(t+2, q)]))
            for i in range(len(self.groups[(t, q)])):
                if print_progress:
                    print(repr(i+1) + "/" + repr(len(self.groups[(t, q)])))
                (s, x) = self.groups[(t, q)][i]
                for (m1, m2, f) in vertices_2out(s):
                    (typ1, a1, b1, c1, d1, bij1) = self.edges[(s, m1)]
                    (typ2, a2, b2, c2, d2, bij2) = self.edges[(m1, f)]
                    if typ1 == 's' and typ2 == 'm':
                        if (c1, d1) == (a2, b2):
                            if x[a1] == 0:
                                val = list(len(x) * [0])
                                for dummy in list(range(0, a1)) + list(range(a1+1, len(x))):
                                    val[bij2[bij1[dummy]]] = x[dummy]
                                val[c2] = 1
                                j = self.groups[(t+2, q)].index((f, tuple(val)))
                                lb[i, j] = 1
            self.lb_stored[(t, q)] = lb

        # Left-multiply: u @ lb
        u_arr = np.asarray(u, dtype=np.int32)
        return (u_arr @ lb._m.astype(np.int32) % 2).astype(np.uint8)

    def chainsq1(self, u, tq, print_progress=True):
        """
        u a cycle in bigrading (t,q); returns numpy uint8 array (cycle in (t+1,q)).

        Uses a 2x integer trick to avoid fractions: we compute 2*l[j] as
        integers (guaranteed even since u is a GF(2) cycle), then divide by 2
        and reduce mod 2.
        """
        (t, q) = tq
        n_next = len(self.groups[(t+1, q)])
        l2 = [0] * n_next   # will hold 2 * l[j]

        if print_progress:
            print("Finding 1-sign assignments.")
        if (t, q) in self.sign1_stored:
            sign1 = self.sign1_stored[(t, q)]
        else:
            sign1 = dict()
            for i in range(len(self.groups[(t, q)])):
                temp = [j for j in range(n_next) if self.maps[(t, q)][i, j] == 1]
                for j in temp:
                    first = first_difference(self.groups[(t, q)][i][0],
                                             self.groups[(t+1, q)][j][0])
                    sign1[(i, j)] = sum(self.groups[(t, q)][i][0][:first]) % 2
            self.sign1_stored[(t, q)] = sign1

        for (i, j) in sign1.keys():
            if u[i] == 1:
                if sign1[(i, j)] == 0:
                    l2[j] += 1
                else:
                    l2[j] -= 1

        # l[j] = l2[j]/2 is guaranteed to be an integer (u is a GF(2) cycle)
        v = np.array([(l2[j] // 2) % 2 for j in range(n_next)], dtype=np.uint8)
        return v

    def chainsq2(self, u, tq, print_progress=True, type='oi'):
        """
        u is cycle in bigrading (t,q); returns numpy uint8 array (cycle in (t+2,q)).
        """
        (t, q) = tq
        v = np.zeros(len(self.groups[(t+2, q)]), dtype=np.uint8)

        if print_progress:
            print("Finding 1-sign assignments.")
        if (t, q) in self.sign1_stored:
            sign1 = self.sign1_stored[(t, q)]
        else:
            sign1 = dict()
            n_next = len(self.groups[(t+1, q)])
            for i in range(len(self.groups[(t, q)])):
                temp = [j for j in range(n_next) if self.maps[(t, q)][i, j] == 1]
                for j in temp:
                    first = first_difference(self.groups[(t, q)][i][0],
                                             self.groups[(t+1, q)][j][0])
                    sign1[(i, j)] = sum(self.groups[(t, q)][i][0][:first]) % 2
            self.sign1_stored[(t, q)] = sign1

        if print_progress:
            print("Finding 2-sign assignments.")
        if (t, q) in self.sign2_stored:
            sign2 = self.sign2_stored[(t, q)]
        else:
            start_time = time.time()
            sign2 = dict()
            for i in range(len(self.groups[(t, q)])):
                if print_progress:
                    print(repr(i+1) + "/" + repr(len(self.groups[(t, q)])))
                vector_1 = self.maps[(t, q)][i]   # i-th row (numpy array)

                (s_state, s_gen) = self.groups[(t, q)][i]
                final_full_list = vertices_2out(s_state)
                final_list = [f for (m1, m2, f) in final_full_list]
                temp_jlist = []

                for j in range(len(self.groups[(t+2, q)])):
                    (f_state, f_gen) = self.groups[(t+2, q)][j]
                    if f_state in final_list:
                        m_state = final_full_list[final_list.index(f_state)][0]
                        (a1, b1, c1, d1, bij1) = self.edges[(s_state, m_state)][1:]
                        (a2, b2, c2, d2, bij2) = self.edges[(m_state, f_state)][1:]
                        temp_key = [s_gen[stuff] == f_gen[bij2[bij1[stuff]]]
                                    for stuff in bij1.keys()
                                    if bij1[stuff] not in [a2, b2]]
                        if False not in temp_key:
                            temp_jlist.append(j)

                for j in temp_jlist:
                    vector_2 = self.maps[(t+1, q)].transpose()[j]   # j-th column
                    list_intermediates = [
                        k for k in range(len(vector_1))
                        if vector_1[k] == 1 and vector_2[k] == 1
                    ]

                    if len(list_intermediates) != 0:
                        temp = dict()
                        sign2[(i, j)] = [0, temp]

                        if len(list_intermediates) == 2:
                            first  = first_difference(self.groups[(t, q)][i][0],
                                                      self.groups[(t+2, q)][j][0])
                            second = second_difference(self.groups[(t, q)][i][0],
                                                       self.groups[(t+2, q)][j][0])
                            sign2[(i, j)][0] = (
                                sum(self.groups[(t, q)][i][0][:first]) *
                                sum(self.groups[(t, q)][i][0][first+1:second])
                            ) % 2
                            sign2[(i, j)][1][list_intermediates[0]] = list_intermediates[1]
                            sign2[(i, j)][1][list_intermediates[1]] = list_intermediates[0]

                        if len(list_intermediates) == 4:
                            [start_gen, mid_gen, final_gen] = [
                                self.groups[(t, q)][i],
                                [self.groups[(t+1, q)][list_intermediates[x]] for x in range(4)],
                                self.groups[(t+2, q)][j]
                            ]
                            [start_state, mid_state, final_state] = [
                                start_gen[0],
                                [mid_gen[x][0] for x in range(4)],
                                final_gen[0]
                            ]
                            [start_onoff, mid_onoff] = [
                                start_gen[1],
                                [mid_gen[x][1] for x in range(4)]
                            ]
                            [start_resolution, mid_resolution] = [
                                self.vertices[start_state],
                                [self.vertices[mid_state[x]] for x in range(4)]
                            ]

                            start_index = self.edges[(start_state, mid_state[0])][1]
                            mid_index = []
                            for x in range(4):
                                guess = self.edges[(start_state, mid_state[x])][3]
                                if mid_onoff[x][guess] == 0:
                                    guess = self.edges[(start_state, mid_state[x])][4]
                                mid_index.append(guess)

                            [start_circle, mid_circle] = [
                                start_resolution.vertex_incidences[start_index],
                                [mid_resolution[x].vertex_incidences[mid_index[x]]
                                 for x in range(4)]
                            ]

                            edges_lb = [
                                first_difference(start_state, final_state) + 1,
                                second_difference(start_state, final_state) + 1
                            ]

                            (matching, correction) = ladybug(start_circle, mid_circle, edges_lb)
                            sign2[(i, j)][0] = correction
                            for x in range(4):
                                sign2[(i, j)][1][list_intermediates[x]] = (
                                    list_intermediates[matching[x]]
                                )

            self.sign2_stored[(t, q)] = sign2
            if print_progress:
                print("Time taken is " + repr(time.time() - start_time))

        if print_progress:
            print("Making the local matching choice.")
        choice = [dict() for x in range(len(self.groups[(t+1, q)]))]
        for x in range(len(self.groups[(t+1, q)])):
            local_state = [0, 0]
            for i in range(len(self.groups[(t, q)])):
                if u[i] == 1 and self.maps[(t, q)][i, x] == 1:
                    if local_state[0] == 0:
                        local_state[0] = 1
                        local_state[1] = i
                    else:
                        local_state[0] = 0
                        choice[x][i] = (0, local_state[1])
                        if sign1[(i, x)] != sign1[(local_state[1], x)]:
                            choice[x][local_state[1]] = (0, i)
                        else:
                            choice[x][local_state[1]] = (1, i)

        if print_progress:
            print("Finding Sq^2 at the chain level.")
        ssign2 = deepcopy(sign2)
        for (i, j) in list(ssign2.keys()):
            if u[i] == 0:
                del ssign2[(i, j)]

        while len(ssign2) != 0:
            (i, j) = next(iter(ssign2.keys()))
            inter = next(iter(ssign2[(i, j)][1].keys()))

            v[j] = (v[j] + 1) % 2

            while doubletest_keys((i, j), inter, ssign2, 1):
                if len(ssign2[(i, j)][1]) == 2:
                    v[j] = (v[j] + ssign2[(i, j)][0]) % 2
                if type == 'io':
                    if len(ssign2[(i, j)][1]) == 4:
                        v[j] = (v[j] + 1) % 2
                inter_next = ssign2[(i, j)][1][inter]
                del ssign2[(i, j)][1][inter]
                del ssign2[(i, j)][1][inter_next]
                inter = inter_next
                if len(ssign2[(i, j)][1]) == 0:
                    del ssign2[(i, j)]

                v[j] = (v[j] + choice[inter][i][0]) % 2
                i = choice[inter][i][1]

        return v

###############################################################################


############# Some other (useful) functions. ##################################

def print_rank(answer):
    """Pretty-print ranks of a bigraded homology theory.

    Input: list of ((t,q), rank) tuples.
    """
    if not answer:
        print("(empty)")
        return
    t_min = min(x[0][0] for x in answer)
    t_max = max(x[0][0] for x in answer)
    q_min = min(x[0][1] for x in answer)
    q_max = max(x[0][1] for x in answer)

    nrows = (q_max - q_min) // 2 + 2
    ncols = t_max - t_min + 2

    M = [[0] * ncols for _ in range(nrows)]
    for i in range(1, ncols):
        M[0][i] = i + t_min - 1
    for i in range(1, nrows):
        M[i][0] = -2*i + q_max + 2
    for ((t, q), r) in answer:
        M[(q_max - q)//2 + 1][t - t_min + 1] = r

    # Format as a grid
    col_widths = [max(len(str(M[r][c])) for r in range(nrows)) for c in range(ncols)]
    lines = []
    for r in range(nrows):
        parts = []
        for c in range(ncols):
            val = str(M[r][c])
            if r > 0 and c > 0 and M[r][c] == 0:
                val = '.'
            parts.append(val.rjust(col_widths[c]))
        lines.append(' '.join(parts))
    print('\n'.join(lines))


def define_knot(name, print_progress=True):
    import builtins
    frame = __import__('sys')._getframe(1)
    frame.f_globals[name] = KnotKh(
        eval(name + '_resolved_knot', frame.f_globals),
        eval(name + '_cross', frame.f_globals),
        print_progress
    )


def compute_knot(name, print_progress=False):
    import sys
    frame = sys._getframe(1)
    frame.f_globals['Central_Knot'] = KnotKh(
        frame.f_globals[name + '_resolved_knot'],
        frame.f_globals[name + '_cross'],
        print_progress
    )

"""
DOES NOT WORK WITH UNKNOT OR HOPF LINK COMPONENTS
"""
