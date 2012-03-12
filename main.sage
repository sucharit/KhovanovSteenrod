## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Contact: lipshitz@math.columbia.edu, sucharit@math.columbia.edu

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


#print "Loading KhovanovSteenrod.sage version 2.0"

from sage.plot.colors import rainbow
import time
from sage.interfaces.chomp import homchain
#from UserTuple import UserTuple


################## Some easy (but convenient) functions. ######################

def first_difference(a,b):
    "Returns index of first place tuples a and b differ."
    if a[0] != b[0]:
        return 0
    return first_difference(a[1:],b[1:])+1

def second_difference(a,b):
    "Returns index of second place tuples a and b differ."
    first = first_difference(a,b)
    return first_difference(a[first+1:],b[first+1:])+first+1

def second_index(data, elt):
    "Returns the second place that elt occurs in data"
    return data[data.index(elt)+1:].index(elt)+data.index(elt)+1

def doubletest_keys(a,b,dic,p=0):
    """dic is a dictionary. If p=0, checks if a is in dic, and b is in dic[a] (if p=1, checks if a in dic and b in dic[a][1])
    (the second test only makes sense if the first passes, hence this precaution)"""
    if a in dic.keys():
        if p==0 and b in dic[a]:
            return True
        if p==1 and b in dic[a][1]:
            return True
    return False

def occurs_twice(data, elt):
    #Checks if +-elt occurs twice in a list.
    if len([x for x in data if (x == elt or x == -elt)])>1:
        return True
    return False
    
def find_cycle_rep(u,V,W):
    """u is a vector in W, which is a quotient space of V. returns a vector v in V, which is a lift of u."""
    z=W.zero()#needed for quotienting

    #find the quotient map matrix (matrices act on row vectors - Sage convention!).
    M=matrix(GF(2),V.rank(),W.rank(),sparse=True) #RL: isn't this built in somewhere?
    for i in range(V.rank()):
        M[i]=W.coordinate_vector(V.basis()[i]+z)

    solution=M.solve_left(W.coordinate_vector(u))#returns coordinates of a solution
    return V.linear_combination_of_basis(solution)

def hypercube(n,faces=False):
    """
    Inputs: n - integer - dimension of the hypercube
    Returns: triple (V,E,F) where V is a list of vertices of the hypercube and E is a list of edges -- pairs of elements of V (initial,final).
    and F is a list of 2D faces -- 4-tuples of elements of V (initial, mid1, mid2, final) (if faces=True)
    """
    if n==1 and faces ==True:
        return ([(0,),(1,)],[((0,),(1,))],[])
    if n==1 and faces ==False:
        return ([(0,),(1,)],[((0,),(1,))])
    if n==2 and faces ==True:
        return ([(0, 0), (1, 0), (0, 1), (1, 1)], [((0, 0), (1, 0)), ((0, 1), (1, 1)), ((0, 0), (0, 1)), ((1, 0), (1, 1))],[((0, 0), (1, 0), (0, 1), (1, 1)),])
    if faces==True:
        (nMinusOneVerts,nMinusOneEdges,nMinusOneFaces) = hypercube(n-1,True)
    else:
        (nMinusOneVerts,nMinusOneEdges) = hypercube(n-1,False)
    vertices = [vert+(0,) for vert in nMinusOneVerts]+[vert+(1,) for vert in nMinusOneVerts]
    edges = [(edge[0]+(0,),edge[1]+(0,)) for edge in nMinusOneEdges]+[(edge[0]+(1,),edge[1]+(1,)) for edge in nMinusOneEdges]+[(vertex+(0,),vertex+(1,)) for vertex in nMinusOneVerts]
    if faces==False:
        return (vertices, edges)
    else:
        faces = [(face[0]+(0,),face[1]+(0,),face[2]+(0,),face[3]+(0,)) for face in nMinusOneFaces]+[(face[0]+(1,),face[1]+(1,),face[2]+(1,),face[3]+(1,)) for face in nMinusOneFaces]+[(edge[0]+(0,),edge[1]+(0,),edge[0]+(1,),edge[1]+(1,)) for edge in nMinusOneEdges]
        return (vertices, edges, faces)

def vertices_out(v):
    "Returns a list of the vertices at the ends of outgoing edges from the vertex v."
    answer = [v[0:i]+(1,)+v[i+1:] for i in range(len(v)) if v[i]==0]
    return answer

def vertices_2out(v):
    "Returns a list of triples (m1,m2,f) where f is at the end of an outgoing 2D face from the vertex v, and m1 and m2 are the two intermediate vertices."
    ans=[]
    for i in range(len(v)):
        if v[i]==0:
            m1=v[0:i]+(1,)+v[i+1:]
            for j in range(i+1,len(v)):
                if v[j]==0:
                    m2=v[0:j]+(1,)+v[j+1:]
                    f=m1[0:j]+(1,)+m1[j+1:]
                    ans=ans+[(m1,m2,f),]
    return ans

###############################################################################


############## The class CycList and some related functions. ##################

class CycList(): #Should really make this extend tuple.
    """
    Input:
    data: a list or tuple of integers.

    Methods:
    rotate(self, n): returns the cyclic list gotten by rotating self n units.
    reflect(self): returns the cyclic list gotten by reflecting self and negating all the entries.
    merge and split: merging two cylists, or splitting one.
    """

    def __init__(self, data):
        self.data = data
        #UserTuple.__init__(data)
    
    def __hash__(self):
        return hash(self.data)

    def __eq__(self, other):
        return self.data == other.data

    def __repr__(self):
        return "CL("+repr(self.data)+")"

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[item]

    def rotate(self, n):
        "Return result of rotating leftwards by n units."
        return CycList(self.data[n:]+self.data[0:n])

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
            #The sorted things could be equal without cyclists being equal. viz. the (2,n) torus links
            md = list(self.data)
            if nmd == od:
                md = list(self.reflect().data)
            od = list(other.data)
            #Modulo cyclic permutations, our md is now either same as od or same as reverse of od.
            if od[(od.index(md[0])+1)%(len(od))] == md[1%(len(od))]:#Just checks two consecutive elements.
                return True
            else:
                return False
            #This test still fails if we have an unknot or a Hopf link, or if a CycList has 1 element.

    def merge(self,other,n):
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
        answer = new_a.data[a_loc+1:]+new_a.data[0:a_loc]+[-n,]+new_b.data[b_loc+1:]+new_b.data[0:b_loc]+[-n,]
        return CycList(answer)

    def split(self,n):
        "Splits cyclic list self along the edge labeled n. Returns (b,c) of CycList's"
        adata = list(self.data)
        newn = n
        if not (n in adata):
            newn = -n
        first = adata.index(newn)
        second = first + adata[first+1:].index(newn) + 1
        b = CycList(adata[first+1:second]+[-newn,])
        c = CycList(adata[second+1:]+adata[0:first]+[-newn,])
        return (b,c)


def equiv(a,b):
    return a.equiv(b)

def merge(a,b,n):
    return a.merge(b,n)

def split(a,n):
    return a.split(n)


def ladybug(start, mid, edges):
    """
    start is a cyclist, mid is a list of 4 cyclists, and edges is a list of two edges, whose endpoints are linked S^0-s on start.
    the 4 cyclists in mid are the 4 circles obtained by splitting start along the 2 edges
    start is in a ladybug configuration, so we want to output the matching on mid
    output is a pair ((self-matching on 0,1,2,3), 0 or 1) where 0 if matching was the right one, 1 if it was the wrong one.
    """

    answer=dict()

    a = abs(edges[0])
    b = abs(edges[1])


    #Rotate / reflect so that start looks like (a,?,-b,?,a,?,-b,?)
    if a in start.data:
        tempns = CycList(start.data)#extra caution
    else:
        tempns = start.reflect()
        
    ns = tempns.rotate(tempns.data.index(a)) #ns = "new start"


    #Adjust the mids so that a and b (which are positive) both occur in all mids, and a is the first entry.
    nmids = list()
    for x in mid:
        if a in x.data:
            nmids.append(x.rotate(x.data.index(a)))
        else:
            temp = x.reflect()
            nmids.append(temp.rotate(temp.data.index(a)))            


    #Case that one of the positively oriented intervals is nonempty.
    if ns.data[1] != -b or ns.data[second_index(ns.data,a)+1] != -b:
        #Rotate some more, so that in (a,?_1,-b,?,a,?,-b,?), ?_1 is non-empty.
        if ns.data[1] == -b:
            ns = ns.rotate(second_index(ns.data,a))
        c = ns.data[1]
        #Now pair up the elts of mid
        matched_pair = [x for x in nmids if ( (x.data[1]==c) or (x.data[len(x.data)-1]==-c) )]
        answer[nmids.index(matched_pair[0])]=nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])]=nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if not (i in answer.keys())]
        answer[unmatched[0]]=unmatched[1]
        answer[unmatched[1]]=unmatched[0]
        return (answer, 0)

    #Case that both positively oriented intervals are empty.    
    else:
        #Rotate so that in (a,-b,?_1,a,-b,?_2), ?_1 is non-empty.
        if second_index(ns.data,a) == 2:
            ns = ns.rotate(second_index(ns.data,a)) #i.e., by 2...
        c = ns.data[2]
        matched_pair = [x for x in nmids if (x.data[(x.data.index(b)+1)%(len(x.data))]==c) or (x.data[len(x.data)-2]==-c)]#extra precaution with len
        answer[nmids.index(matched_pair[0])]=nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])]=nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if not (i in answer.keys())]
        answer[unmatched[0]]=unmatched[1]
        answer[unmatched[1]]=unmatched[0]
        return (answer,1)

###############################################################################


############### The class ResolvedKnot and related function. ##################

class ResolvedKnot():
    """
    Inputs:
     --vertex_incidences: a list of CycList's
     --

    Methods:
    --flip(edge): Return resolved knot obtained by switching resolution at edge. edge is a number from 1 to n.
    --graph(return_colors = False): Returns a graph or a pair (graph, edge_colors), representing self. Currently ignores how edges sit in the plane.

    Data:
    N = number of edges (=number of crossings in unresolved knot)
    C = number of components
    """
    #Add some code so you can call it with lists rather than CycList's.
    def __init__(self, vertex_incidences, label = []):
        self.vertex_incidences = vertex_incidences
        self.N = len(sum([x.data for x in vertex_incidences],[]))/2
        if label != []:
            self.label = tuple(label)
        if label == []:
            self.label = self.N*(0,)
        self.edge_incidences = dict([(i+1,list()) for i in range(self.N)])
        for x in vertex_incidences:
            for y in x.data:
                self.edge_incidences[abs(y)].append(vertex_incidences.index(x))
        self.C = len(self.vertex_incidences)

    def __repr__(self):
        return "ResolvedKnot with label: "+repr(self.label)+"\n     v.i.: "+repr(self.vertex_incidences)+"\n"

    def flip(self, edge):
        "Returns the resolved knot obtained by switching the resolution at edge. edge is a number from 1 to n."
        new_label = self.label[0:edge-1]+((self.label[edge-1]+1)%2, )+self.label[edge:]
        if self.edge_incidences[edge][0]==self.edge_incidences[edge][1]: #This is the case of a split.
            sc = self.edge_incidences[edge][0] #This circle is being split.
            other_components = self.vertex_incidences[0:sc]+self.vertex_incidences[sc+1:]
            new_vertex_incidences = tuple(other_components) + split(self.vertex_incidences[sc],edge)
        if self.edge_incidences[edge][0]!=self.edge_incidences[edge][1]: #This is the case of a merge.
            mc1 = min((self.edge_incidences[edge][0],self.edge_incidences[edge][1]))
            mc2 = max((self.edge_incidences[edge][0],self.edge_incidences[edge][1]))
            other_components = self.vertex_incidences[0:mc1]+self.vertex_incidences[mc1+1:mc2]+self.vertex_incidences[mc2+1:]
            new_vertex_incidences = tuple(other_components) + (merge(self.vertex_incidences[mc1],self.vertex_incidences[mc2],edge),)
        return ResolvedKnot(new_vertex_incidences, new_label)

    def graph(self, return_colors = False):
        """Returns a graph representing self. 
        Currently ignores whether edges are supposed to be inside or outside circles.
        If optional value return_colors is True, returns a pair (graph, edge_colors). 
        Then, to display, use: 
        sage: (graph, pretty_colors) = my_knot.graph(return_colors = True)
        sage: graph.plot(edge_colors = pretty_colors)
        The same information is encoded in:
        sage: my_knot.graph().plot(edge_labels=True)

        Note: vertex labels are triples (name of this vertex, number of cycle vertex lives on, index of this vertex on that cycle).
        """
        G = Graph(multiedges = True)
        R = rainbow(2)
        abs_vertex_incidences = [[abs(i) for i in x.data] for x in self.vertex_incidences]
        edge_colors = dict()
        edge_colors[R[0]]=list()
        edge_colors[R[1]]=list()
        for x in abs_vertex_incidences:
            for j in range(len(x)):
                G.add_vertex(pvert((x[j],abs_vertex_incidences.index(x),j)))
                new_edge = (pvert((x[j],abs_vertex_incidences.index(x),j)),pvert((x[(j+1)%len(x)],abs_vertex_incidences.index(x),(j+1)%len(x))),None)
                G.add_edge(new_edge)
        for i in range(self.N):
            if self.edge_incidences[i+1][0] != self.edge_incidences[i+1][1]:
                new_edge = (pvert((i+1,self.edge_incidences[i+1][0], abs_vertex_incidences[self.edge_incidences[i+1][0]].index(i+1))),pvert((i+1,self.edge_incidences[i+1][1],abs_vertex_incidences[self.edge_incidences[i+1][1]].index(i+1))))
            else:
                new_edge = (pvert((i+1,self.edge_incidences[i+1][0],abs_vertex_incidences[self.edge_incidences[i+1][0]].index(i+1))),pvert((i+1,self.edge_incidences[i+1][1],abs_vertex_incidences[self.edge_incidences[i+1][1]][self.edge_incidences[i+1][0].index(i+1):].index(i+1))))
            G.add_edge(new_edge,label=i+1)
        edge_colors[R[0]] = [x for x in G.edges() if x[2]==None]
        edge_colors[R[1]] = [x for x in G.edges() if x[2]!=None]
        if return_colors:
            return (G,edge_colors)
        return G

class pvert():
    "A class to help me make the vertex labels in ResolvedKnot.graph() pretty."
    def __init__(self, data):
        self.data = data
        
    def __hash__(self):
        return hash(self.data)

    def __repr__(self):
        return repr(self.data[0])

    def __eq__(self, other):
        return self.data == other.data

    def __neq__(self, other):
        return not self.__eq__(other)


###############################################################################




############ Now some main functions. #########################################


def populate_vertices(zero_res):
    "Returns: dictionary with keys elements of hypercube(N)[0] and values the corresponding ResolvedKnot"
    N = zero_res.N
    cube = dict()
    cube[N*(0,)]=zero_res
    for s in hypercube(N)[0][1:]:
        i = list(s).index(1)
        prev_s = s[0:i]+(0,)+s[i+1:]
        cube[s]=cube[prev_s].flip(i+1)
    return cube

def find_fly_in_soup(fly,soup):
    "Returns the index of the CycList fly in the list of CycList's soup"
    for b in soup:
        if equiv(fly,b):
            return list(soup).index(b)
    raise(Exception("Can't find fly "+repr(fly)+" in soup "+repr(soup)))

def populate_edges(cube):
    """
    Input: a dictionary with keys vertices of N-D hypercube (i.e., hypercube(N)[0]) and values the corresponding ResolvedKnot. i.e., the output of populate_vertices().
    Returns: dictionary with keys edges of N-D hypercube (i.e., hypercube(N)[1]) and values: 
    --('m' or 's', a, b, c, d, bijection) where a,b,c,d are the indices of the three exceptional components (obvious one repeated),
    ----'m' or 's' says whether it's a merge or a split.
    ----bijection is a bijection between the remaining (non-exceptional) components, encoded as a dictionary.
    Probably fails in some degenerate cases: unknot, Hopf link components.
    """
    #Loop through the edges:
    N = len(cube.keys()[0])
    answer = dict()
    for e in hypercube(N)[1]:
        #merge or split?
        initial_vert = cube[e[0]]
        final_vert = cube[e[1]]
        edge_flipping = first_difference(e[0],e[1])+1 #the knot edge being changed
        if initial_vert.C>final_vert.C:
            move_type = 'm'
        else:
            move_type = 's'
        #Find the exceptional components:
        a=min(initial_vert.edge_incidences[edge_flipping][0],initial_vert.edge_incidences[edge_flipping][1])
        b=max(initial_vert.edge_incidences[edge_flipping][0],initial_vert.edge_incidences[edge_flipping][1])
        c=min(final_vert.edge_incidences[edge_flipping][0],final_vert.edge_incidences[edge_flipping][1])
        d=max(final_vert.edge_incidences[edge_flipping][0],final_vert.edge_incidences[edge_flipping][1])

        #Find the bijection between the non-exceptional components:
        bijection = dict()
        #This part of the code is pretty inefficient.
        for x in range(initial_vert.C):
            if ( x != a and x != b ):
                bijection[x] = find_fly_in_soup(initial_vert.vertex_incidences[x],final_vert.vertex_incidences)
        answer[e]=(move_type, a, b, c, d, bijection)
    return answer

def chain_groups(cube,(t_correct,q_correct)=(0,0)):
    """
    Input: the output of populate_vertices.
    Returns: A dictionary whose keys are pairs of integers (the bi-grading), and the value is a list of all generators in that bigrading.
    """
    generators = dict()
    for s in cube.keys():
        for g in hypercube(cube[s].C)[0]:
            t = sum(s)+t_correct
            q = len(g)-2*sum(g)+t+q_correct-t_correct
            if (t,q) in generators.keys():
                generators[(t,q)].append((s,g))
            else:
                generators[(t,q)]=[(s,g),]
    return generators

def chain_maps(group,edges, print_progress = True):
    """
    Input: the output of chain_group and the output of populate_edges
    Output: A dictionary whose keys are pairs of integers (t,q), and the key is a matrix which describes the map from the chain group at (t,q) to the chain group at (t+1,q).
    This might be the slowest step!!
    Matrix acts like: input row vector times matrix = output row vector. (v\cdot M) (Sage convention!!)
    """
    maps = dict()
    #Populate matrices of correct dimensions with zeroes.
    for k in group.keys():
        if print_progress:
            print repr(1+group.keys().index(k))+'/'+repr(len(group.keys()))
        l=(k[0]+1,k[1])#Added l to make stuff more human-readable. SS
        if l in group.keys():
            maps[k]=matrix(GF(2),len(group[k]),len(group[l]),sparse=True)
    #Now, add 1's for corresponding to the edge maps...
    #RL just wrote this, and has done zero testing so far. Promises to debug later.
        for i in range(len(group[k])): #s is a state; x is a generator of Kh of resolution s.
            (s,x) = group[k][i]
            for t in vertices_out(s): #states pointed to by s
                (edge_type, a, b, c, d, bijection) = edges[(s,t)]
                if edge_type == 'm':
                    if x[a]!=1 or x[b]!=1:
                        val = (len(x)-1)*[0,] #Val will be f_{this edge}
                        for dummy in range(0,a)+range(a+1,b)+range(b+1,len(x)):#Needed new variable name.
                            val[bijection[dummy]] = x[dummy]
                        val[c] = x[a]+x[b]
                        j = group[l].index((t,tuple(val)))
                        maps[k][i,j] = 1 # No need to add 1 to the existing value, since Kh complex is nice that way. 
                if edge_type == 's':
                    val = (len(x)+1)*[0,] #Val will be f_{this edge}
                    for dummy in range(0,a)+range(a+1,len(x)):#Again new variable name.
                        val[bijection[dummy]] = x[dummy]
                    if x[a]==1:
                        val[c] = 1
                        val[d] = 1
                        j = group[l].index((t,tuple(val)))
                        maps[k][i,j] = 1 
                    if x[a]==0:#Two hits
                        val[c] = 0
                        val[d] = 1
                        j = group[l].index((t,tuple(val)))
                        maps[k][i,j] = 1 
                        val[c] = 1
                        val[d] = 0
                        j = group[l].index((t,tuple(val)))
                        maps[k][i,j] = 1
    return maps

def complexes(groups,maps, print_progress = True):
    """
    Input: Output of chain_groups and chain_maps.
    (Don't really need chain_groups. Just need dimensions (no way to recover from maps).
    Outputs dictionary whose keys are the keys for q-gradings, and value is chain complex.
    """
    answer=dict()
    homans=dict()

    temp=dict()
    for (t,q) in groups.keys():
        if q in temp.keys():
            temp[q][0]=min(temp[q][0],t)
            temp[q][1]=max(temp[q][1],t)
        else:
            temp[q]=[t,t]
    
    
    for i in range(len(temp.keys())):
        if print_progress:
            print repr(i+1)+"/"+repr(len(temp.keys()))
        q=temp.keys()[i]
        ans_dict=dict()
        for t in range(temp[q][0],temp[q][1]+1):
            if (t,q) in groups.keys() and (t+1,q) in groups.keys():

                ans_dict[t]=maps[(t,q)].transpose() #Sage reverses their convention while defining chain complexes.

            elif (t,q) in groups.keys():
                ans_dict[t]=matrix(GF(2),0,len(groups[(t,q)]),sparse=True)
            elif (t+1,q) in groups.keys():
                ans_dict[t]=matrix(GF(2),len(groups[(t+1,q)]),0,sparse=True)
            else:
                ans_dict[t]=matrix(GF(2),0,0,sparse=True)

        answer[q]=ChainComplex(ans_dict,base_ring=GF(2),check_products=False)

        start_time=time.time()
        homans[q]=homchain(answer[q],generators=True)
        for t in homans[q].keys():#formatting issue. if rank=0, homchain only returns a vector space.
            if type(homans[q][t])!=tuple:
                homans[q][t]=(homans[q][t],[])
        if print_progress:
            print repr(q)+"::: time "+repr(time.time()-start_time)


    return (answer,homans)


def homology(groups,maps, print_progress = True):
    """
    Input: Output of chain_groups and chain_maps.
    (Don't really need chain_groups. Just need dimensions (no way to recover from maps).
    Outputs dictionary whose keys are the keys for groups, and value is homology information.
    """
    answer=dict()
    for (t,q) in groups.keys():
        if print_progress:
            print repr(1+groups.keys().index((t,q)))+"/"+repr(len(groups.keys()))
        #compute kernel
        if (t,q) in maps.keys():
            ker = kernel(maps[(t,q)])
        else:
            dummy = matrix(GF(2),len(groups[(t,q)]),1,sparse=True)
            ker = kernel(dummy)
        #compute image
        if (t-1,q) in maps.keys():
            im = maps[(t-1,q)].row_space()
        else:
            dummy = matrix(GF(2),1,len(groups[(t,q)]),sparse=True)
            im = dummy.row_space()
        #compute quotient
        try:
            ker/im
            answer[(t,q)]=(ker,im,ker/im)#kernel, image and quotient, as vector spaces.
        except ArithmeticError:
            """
            This case happens with (t,q)=(-1,-1) in R8_18 (the presentation from knotsupto10.sage).
            However, d^2 is still zero. It turns out that the implementation of kernel in Sage (<=4.7.1) is wrong!!!!
            """
            print "Something is seriously wrong"
            print (t,q) 

    return answer


###############################################################################



############### The main KnotKh class #########################################

class KnotKh():
    """
    Input:
    ResolvedKnot representing the 0-resolution

    Data:
    --data: copy of resolved knot it's called by.
    --vertices: the elements at the vertices of the Khovanov cube. (result of call to populate_vertices())
    --edges: maps for the edges. (result of call to populate_edges())
    --groups: chain groups for Khovanov cube. (result of call to chain_groups())
    --maps: chain maps for Khovanov cube. (result of call to chain_maps())
    --complexes: the chain complexes in different q-gradings. (result of call to complexes())
    --homology: Khovanov homology (result of call to homology())
    --sign1_stored: as much of the 1-sign assignment as we've computed so far.
    --sign2_stored: as much of the 2-sign assignment as we've computed so far.
    --lb_stored: as much of the ladybug map as we've computed so far.

    Methods:
    --d2(): test whether d^2=0 on self.
    --hom_rank(): Return list of ranks of Khovanov homology of self, as tuples ((t,q),r) where (t,q) is the bigrading and r is the rank in that bigrading.
    --print_hom_rank(): prints the results of hom_rank(), nicely formatted. Uses funciton print_rank().
    --sq1((t,q),n=1,print_progress=True): computes Sq1 (the Bockstein)
    --sq2((t,q),n=1,print_progress=True,type='oi'): computes Sq2 for n'th generator in bigrading (t,q) of type oi (out-in) or io (in-out)
    --chainsq2: used by sq2 for the chain-level computations.
    """

    def __init__(self, data, (nplus,nminus)=(0,0), print_progress = True):
        self.data = ResolvedKnot(data.vertex_incidences)#Have to do this strange business. Resolved knot needs _call_ stuff?
        if print_progress:
            print "Popuating vertices."
        self.vertices = populate_vertices(self.data)
        if print_progress:
            print "Populating edges."
        self.edges = populate_edges(self.vertices)
        if print_progress:
            print "Finding chain groups."
        self.groups = chain_groups(self.vertices,(-nminus,nplus-2*nminus))

        start_time=time.time()
        if print_progress:
            print "Finding chain maps."
        self.maps = chain_maps(self.groups,self.edges, print_progress)
        if print_progress:
            print "Time taken is "+repr(time.time()-start_time)

        start_time=time.time()
        if print_progress:
            print "Generating chain complexes."
        (self.complexes,self.newhom)=complexes(self.groups,self.maps,print_progress)
        if print_progress:
            print "Time taken is "+repr(time.time()-start_time)

        
        #start_time=time.time()
        #if print_progress:
        #    print "Computing old homology."
        #self.homology = homology(self.groups,self.maps,print_progress)
        #if print_progress:
        #    print "Time taken is "+repr(time.time()-start_time)

        if print_progress:
            print "Done computing Khovanov homology (over F_2)."

        self.sign1_stored = dict()
        self.sign2_stored = dict()
        self.lb_stored = dict()

        self.hom = dict()
        self.hom_support = []

        self.Sq1 = dict()
        self.Sq2 = dict()


    def d2(self, print_output = False):#checks d2=0
        "Prints the ranks of the products in d2"
        are_we_good = True
        for k in self.maps.keys():
            l=(k[0]+1,k[1])
            if l in self.maps.keys():
                d2=(self.maps[k])*(self.maps[l])
                if d2.rank() !=0:
                    are_we_good = False
                if print_output:
                    print d2.rank()
        return are_we_good

    def hom_rank(self):
        "Returns list of ranks of Khovanov homology of self, as tuples ((t,q),r) where (t,q) is the bigrading and r is the rank in that bigrading."
        ranks = list()
        for q in self.newhom.keys():
            for t in self.newhom[q].keys():
                r = len(self.newhom[q][t][1])
                if r>0:
                    ranks.append(((t,q),r))
        ranks.sort(key=lambda x: (x[0][1],x[0][0]) )
        self.hom = dict(ranks)
        self.hom_support = [(t,q) for ((t,q),r) in ranks]
        return ranks
        
    def print_hom_rank(self):
        "Prints the results of hom_rank(), nicely formatted."
        answer = self.hom_rank()
        print_rank(answer)
        

    def __repr__(self):
        return "KNOT("+repr(self.data)+")"

    def lady_map(self,(t,q),n=1, print_progress = True): 
        """
        Computes ladybug map on bigrading (t,q). Returns True or False, depending on whether the image is zero or not in homology.
        Chooses some basis for homology at (t,q), and computes ladybug on the n-th basis element.
        Returns: True or False.
        """
        if n < 1:
            raise(Exception("Homology groups indexed by positive integers; entry n="+repr(n)+" has no meaning."))
        if n > len(self.newhom[q][t][1]):
            raise(Exception("Our homology group has rank "+repr(len(self.newhom[q][t][1]))+" so can't find ladybug of element "+repr(n)))
        else: #RL: now that the "if" raises an exception, don't need the else here.
            u=self.newhom[q][t][1][n-1] #This is the cycle representative of the homology element. We want to find its lb.
            
            if (t,q) in self.maps.keys() and (t+1,q) in self.maps.keys() and len(self.newhom[q][t+2][1])>0:
                v=self.chainlady_map(u,(t,q), print_progress) #the lady_map computation at chain level. v is a cycle in bigrading (t+2,q).
                return v in self.maps[(t+1,q)].row_space()
            else:
                return True


    def sq1(self,(t,q),n=0, print_progress = True, show_chain=False): 
        """
        Computes Sq1 map on bigrading (t,q). Returns an element in bigrading (t+1,q).
        Chooses some basis for homology at (t,q), and computes Sq1 on the n-th basis element.
        If n=0, computes all.
        """
        if n < 0:
            raise(Exception("Homology groups indexed by positive integers; entry n="+repr(n)+" has no meaning."))
        if n > len(self.newhom[q][t][1]):
            raise(Exception("Our homology group has rank "+repr(len(self.newhom[q][t][1]))+" so can't find Sq1 of element "+repr(n)))
        if n == 0:
            if (t,q) in self.Sq1.keys():
                return self.Sq1[(t,q)]
            else:
                r = len(self.newhom[q][t][1])
                l=[]
                for i in range(1,r+1):
                    l=l+[self.sq1((t,q),i,print_progress=print_progress),]
                self.Sq1[(t,q)]=matrix(GF(2),l)
                return self.Sq1[(t,q)]
        else: 
            u=self.newhom[q][t][1][n-1] #This is the cycle representative of the homology element. We want to find its Sq1.
            
            if (t,q) in self.maps.keys() and len(self.newhom[q][t+1][1])>0:
                v=self.chainsq1(u,(t,q), print_progress=print_progress) #the Sq1 computation at chain level. v is a cycle in bigrading (t+1,q).


                temp_mat=matrix(GF(2),self.newhom[q][t+1][1]+self.maps[(t,q)].rows(),sparse=True)
                v_in_hom=temp_mat.solve_left(v)
                v_in_hom=v_in_hom[:len(self.newhom[q][t+1][1])]

                if show_chain:
                    return (v,v_in_hom)
                else:
                    return v_in_hom


            elif show_chain:
                if (t+1,q) in self.groups.keys():
                    return (vector(GF(2),len(self.groups[(t+1,q)]),sparse=True),vector(GF(2),len(self.newhom[q][t+1][1]),sparse=True))
                else:
                    return (vector(GF(2),0,sparse=True),vector(GF(2),0,sparse=True))
            else:
                if (t+1,q) in self.groups.keys():
                    return vector(GF(2),len(self.newhom[q][t+1][1]),sparse=True)
                else:
                    return vector(GF(2),0,sparse=True)


    def sq2_oi(self,(t,q),n=0,print_progress=True):
        return self.sq2_oi((t,q),n,print_progress=print_progress,type='oi') #to be backward-compatible

    def sq2_io(self,(t,q),n=0,print_progress=True):
        return self.sq2_oi((t,q),n,print_progress=print_progress,type='io') #to be backward-compatible

    def sq2(self,(t,q),n=0, print_progress = True, type='oi',show_chain=False): #RL: isn't it un-natural to start with n=1 rather than n=0? (i.e., propose replacing n-1 by n later, and this 1 by 0.)
        """
        Computes Sq2 (out-in) on bigrading (t,q). Returns an element in bigrading (t+2,q)
        Chooses some basis for homology at (t,q), and computes Sq2 on the n-th basis element.
        Returns: pair (chain representing Sq^2 of input, homology class representing Sq^2 of input)
        """
        if n < 0:
            raise(Exception("Homology groups indexed by positive integers; entry n="+repr(n)+" has no meaning."))
        if n > len(self.newhom[q][t][1]):
            raise(Exception("Our homology group has rank "+repr(len(self.newhom[q][t][1]))+" so can't find Sq2 of element "+repr(n)))
        if n == 0:
            if (t,q) in self.Sq2.keys():
                return self.Sq2[(t,q)]
            else:
                r = len(self.newhom[q][t][1])
                l=[]
                for i in range(1,r+1):
                    l=l+[self.sq2((t,q),i,print_progress=print_progress),]
                self.Sq2[(t,q)]=matrix(GF(2),l)
                return self.Sq2[(t,q)]
        else:
            u=self.newhom[q][t][1][n-1] #This is the chain rep of the homology element. We want to find its Sq2

            if (t,q) in self.maps.keys() and (t+1,q) in self.maps.keys() and len(self.newhom[q][t+2][1])>0:
                v=self.chainsq2(u,(t,q), print_progress=print_progress, type=type) #the Sq2 (out-in) computation at chain level. v is a cycle in bigrading (t+2,q).


                temp_mat=matrix(GF(2),self.newhom[q][t+2][1]+self.maps[(t+1,q)].rows(),sparse=True)
                v_in_hom=temp_mat.solve_left(v)
                v_in_hom=v_in_hom[:len(self.newhom[q][t+2][1])]

                if show_chain:
                    return (v,v_in_hom)
                else:
                    return v_in_hom


            elif show_chain:
                if (t+2,q) in self.groups.keys():
                    return (vector(GF(2),len(self.groups[(t+2,q)]),sparse=True),vector(GF(2),len(self.newhom[q][t+2][1]),sparse=True))
                else:
                    return (vector(GF(2),0,sparse=True),vector(GF(2),0,sparse=True))
            else:
                if (t+2,q) in self.groups.keys():
                    return vector(GF(2),len(self.newhom[q][t+2][1]),sparse=True)
                else:
                    return vector(GF(2),0,sparse=True)


    def chainlady_map(self, u, (t,q), print_progress = True):
        """
        u is cycle in bigrading (t,q)
        returns a cycle representative for ladymap(u) in bigrading (t+2,q)
        makes some matching choice on the way
        """

        if print_progress:
            print "Finding the ladybugs."
        """stores stuff between generators in (t,q) and (t+2,q).
        It is a dictionary, the value at (t,q) being 
        a matrix over GF(2), which acts on generators in grading (t,q) as row vectors."""
        if (t,q) in self.lb_stored.keys():
            lb = self.lb_stored[(t,q)]
        else:
            #Creates new.
            lb = matrix(GF(2),len(self.groups[(t,q)]),len(self.groups[(t+2,q)]),sparse=True)
            for i in range(len(self.groups[(t,q)])):
                if print_progress:
                    print repr(i+1)+"/"+repr(len(self.groups[(t,q)]))
                (s,x)=self.groups[(t,q)][i]
                for (m1,m2,f) in vertices_2out(s):
                    (typ1,a1,b1,c1,d1,bij1)=self.edges[(s,m1)]
                    (typ2,a2,b2,c2,d2,bij2)=self.edges[(m1,f)]
                    if typ1=='s':
                        if typ2=='m':
                            if (c1,d1)==(a2,b2):
                                if x[a1]==0:#Found a ladybug configuration.
                                    val = (len(x))*[0,] #Val will be lb_{this edge}
                                    for dummy in range(0,a1)+range(a1+1,len(x)):#Needed new variable name.
                                        val[bij2[bij1[dummy]]] = x[dummy]
                                    val[c2] = 1
                                    j = self.groups[(t+2,q)].index((f,tuple(val)))

                                    lb[i,j]=1

                            
            self.lb_stored[(t,q)]=lb
                    
        return u*lb

    def chainsq1(self, u, (t,q), print_progress=True):
        """
        u a cycle in bigrading (t,q)
        return a cycle representative for sq1(u) in bigrading (t+1,q)
        """
        l=(len(self.groups[(t+1,q)]))*[0,]

        if print_progress:
            print "Finding 1-sign assignments."
        if (t,q) in self.sign1_stored.keys():
            sign1 = self.sign1_stored[(t,q)]
        else:
            sign1 = dict()
            """stores 1-sign assignments between (t,q) and (t+1,q) as a dictionary.
            the keys are pairs (a,b) where a is a generator in (t,q) and b is a generator in (t+1,q) and there is a flowline from a to b.
            the value is the additive sign, i.e. 0 or 1. NOTE: Additive."""
            for i in range(len(self.groups[(t,q)])):
                temp = [j for j in range(len(self.groups[(t+1,q)])) if self.maps[(t,q)][i,j] == 1]
                for j in temp:
                    first = first_difference(self.groups[(t,q)][i][0],self.groups[(t+1,q)][j][0])
                    sign1[(i,j)]=sum(self.groups[(t,q)][i][0][:first])%2

            self.sign1_stored[(t,q)]=sign1

        ## for j in range(len(self.groups[(t+1,q)])):
        ##     temp_ilist=[i for i in range(len(self.groups[(t,q)])) if (u[i]==1 and (i,j) in sign1.keys())]
        ##     positives=[i for i in temp_ilist if sign1[(i,j)]==0]
        ##     negatives=[i for i in temp_ilist if sign1[(i,j)]==1]
        ##     l[j]=(len(positives)-len(negatives))/2

        for (i,j) in sign1.keys():
            if u[i]==1:
                if sign1[(i,j)]==0:
                    l[j]=l[j]+1/2
                else:
                    l[j]=l[j]-1/2

        v=vector(GF(2),l,sparse=True)

        return v

    def chainsq2(self, u, (t,q), print_progress = True,type='oi'):
        """
        u is cycle in bigrading (t,q)
        returns a cycle representative for sq2(u) in bigrading (t+2,q)
        makes some matching choice on the way
        """
        v=vector(GF(2),(len(self.groups[(t+2,q)]))*[0,],sparse=True)

        
        if print_progress:
            print "Finding 1-sign assignments."
        if (t,q) in self.sign1_stored.keys():
            sign1 = self.sign1_stored[(t,q)]
        else:
            sign1 = dict()
            """stores 1-sign assignments between (t,q) and (t+1,q) as a dictionary.
            the keys are pairs (a,b) where a is a generator in (t,q) and b is a generator in (t+1,q) and there is a flowline from a to b.
            the value is the additive sign, i.e. 0 or 1. NOTE: Additive."""
            for i in range(len(self.groups[(t,q)])):
                temp = [j for j in range(len(self.groups[(t+1,q)])) if self.maps[(t,q)][i,j] == 1]
                for j in temp:
                    first = first_difference(self.groups[(t,q)][i][0],self.groups[(t+1,q)][j][0])
                    sign1[(i,j)]=sum(self.groups[(t,q)][i][0][:first])%2

            self.sign1_stored[(t,q)]=sign1


        if print_progress:
            print "Finding 2-sign assignments."
        """stores stuff between generators in (t,q) and (t+2,q).
        the keys are pairs (a,b), a gen in (t,q), b gen in (t+2,q), such that there are flowlines from a to b.
        the value is a list [s,l], where s is the additive 2-sign assignment (non-zero *only* if there are exactly 2 flowlines)
        and l is a dictionary whose keys are intermediate vertices (2 or 4 keys), and it stores the matching info. (Remember the ladybug rule for oi)."""
        if (t,q) in self.sign2_stored.keys():
            sign2 = self.sign2_stored[(t,q)]
        else:
            #Creates new.
            start_time=time.time()
            sign2 = dict()
            for i in range(len(self.groups[(t,q)])):
                if print_progress:
                    print repr(i+1)+"/"+repr(len(self.groups[(t,q)]))
                vector_1 = self.maps[(t,q)][i]

                (s_state,s_gen)=self.groups[(t,q)][i]
                final_full_list=vertices_2out(s_state)
                final_list=[f for (m1,m2,f) in final_full_list]
                temp_jlist=[]

                for j in range(len(self.groups[(t+2,q)])):
                    (f_state,f_gen)=self.groups[(t+2,q)][j]
                    if f_state in final_list:

                        m_state=final_full_list[final_list.index(f_state)][0]
                        (a1,b1,c1,d1,bij1)=self.edges[(s_state,m_state)][1:]
                        (a2,b2,c2,d2,bij2)=self.edges[(m_state,f_state)][1:]

                        temp_key=[s_gen[stuff]==f_gen[bij2[bij1[stuff]]] for stuff in bij1.keys() if not bij1[stuff] in [a2,b2]]
                        if not False in temp_key:                      
                            temp_jlist.append(j)
                
                for j in temp_jlist:
                    vector_2 = self.maps[(t+1,q)].transpose()[j]
                    intermediates = vector(GF(2),(len(vector_1))*[0,],sparse=True)#This will be the termwise product of vector_1 and vector_2
                    list_intermediates = []#this stores the location of the non-zero terms in intermediates
                    for k in range(len(intermediates)):
                        if vector_1[k] == 1 and vector_2[k] ==1:
                            intermediates[k]=1
                            list_intermediates.append(k)

                    if len(list_intermediates) != 0:#i.e. there are broken flowlines from i to j.
                        temp=dict();sign2[(i,j)]=[0,temp]
                            
                        if len(list_intermediates) == 2:
                            first = first_difference(self.groups[(t,q)][i][0],self.groups[(t+2,q)][j][0])
                            second = second_difference(self.groups[(t,q)][i][0],self.groups[(t+2,q)][j][0])
                            sign2[(i,j)][0]=(sum(self.groups[(t,q)][i][0][:first])*sum(self.groups[(t,q)][i][0][first+1:second]))%2
                            sign2[(i,j)][1][list_intermediates[0]]=list_intermediates[1]
                            sign2[(i,j)][1][list_intermediates[1]]=list_intermediates[0]

                        if len(list_intermediates) == 4:
                            
                            """
                            We are now in ladybug situation. start is the ladybug we start with. mid is a list of 4 elements, which stores the 4 intermediate stuff. final is the final ladybug.
                            """
                            [start_gen,mid_gen,final_gen]=[self.groups[(t,q)][i],[self.groups[(t+1,q)][list_intermediates[x]] for x in range(4)],self.groups[(t+2,q)][j]]#the generators = (s,g)
                            [start_state,mid_state,final_state]=[start_gen[0],[mid_gen[x][0] for x in range(4)],final_gen[0]]#the state s, a vertex of N-d hypercube
                            [start_onoff,mid_onoff]=[start_gen[1],[mid_gen[x][1] for x in range(4)]]#the second part of the generator, for the lack of a better name.
                            [start_resolution,mid_resolution]=[self.vertices[start_state],[self.vertices[mid_state[x]] for x in range(4)]]#the resolved knots
                            
                            start_index=self.edges[(start_state,mid_state[0])][1]#the index of the circle being split
                            mid_index=[]
                            for x in range(4):#the indices of the other relevant circles (we only keep the circle that is "on")
                                guess = self.edges[(start_state,mid_state[x])][3]
                                if mid_onoff[x][guess] == 0:
                                    guess = self.edges[(start_state,mid_state[x])][4]
                                mid_index.append(guess)

                            [start_circle,mid_circle]=[start_resolution.vertex_incidences[start_index],[mid_resolution[x].vertex_incidences[mid_index[x]] for x in range(4)]]
                            """Finally the relevant circles (i.e. cyclists).
                            We will now compute the two edges are are involved to go between i and j"""

                            edges=[first_difference(start_state,final_state)+1,second_difference(start_state,final_state)+1]#Remember, edges are numbered 1 to N.

                            (matching, correction) = ladybug(start_circle, mid_circle, edges)#the output is a dictionary which describes a self-matching of 0,1,2,3
                            sign2[(i,j)][0]=correction
                            for x in range(4):
                                sign2[(i,j)][1][list_intermediates[x]]=list_intermediates[matching[x]]

                            
            self.sign2_stored[(t,q)]=sign2
            if print_progress:
                print "Time taken is "+repr(time.time()-start_time)

        if print_progress:
            print "Making the local matching choice."
        choice = [dict() for x in range(len(self.groups[(t+1,q)]))]#initialisation. Is is really necessary? SS.
        #Stores the local matching. It is a list of dictionaries, where the dictionary entries stores the matchings (incl. *directions*) as a tuple.
        for x in range(len(self.groups[(t+1,q)])):
            local_state=[0,0]
            for i in range(len(self.groups[(t,q)])):
                if u[i] == 1 and self.maps[(t,q)][i,x] == 1:
                    if local_state[0] == 0:
                        local_state[0]=1
                        local_state[1]=i
                    else:
                        local_state[0]=0
                        choice[x][i]=(0,local_state[1])
                        if sign1[(i,x)] != sign1[(local_state[1],x)]:#i.e. the matching is undirected
                            choice[x][local_state[1]]=(0,i)
                        else: #i.e. the matching is directed
                            choice[x][local_state[1]]=(1,i)

        """Now the main definition of v"""

        if print_progress:
            print "Finding Sq^2 at the chain level."
        ssign2=deepcopy(sign2)#we are going to delete keys from sign2, and we are scared.
        for (i,j) in ssign2.keys():
            if u[i] == 0:
                del ssign2[(i,j)]#only interested in things starting at the chain u
        while len(ssign2.keys())!=0:

            (i,j)=ssign2.keys()[0]
            inter=ssign2[(i,j)][1].keys()[0]#this is an intermediate vertex between i and j

            v[j]=v[j]+1#we start a new cycle (in the graph corresponding to j)

            while doubletest_keys((i,j), inter, ssign2,1):

                #v[j]=v[j]+ssign2[(i,j)][0]#the ssign2 contribution
                if len(ssign2[(i,j)][1].keys())==2:
                    v[j]=v[j]+ssign2[(i,j)][0]
                if type=='io':
                    if len(ssign2[(i,j)][1].keys())==4:
                        v[j]=v[j]+1 #since we are doing the "other" ladybug matching                
                inter_next=ssign2[(i,j)][1][inter]#moves according the global matching
                del ssign2[(i,j)][1][inter]
                del ssign2[(i,j)][1][inter_next]
                inter=inter_next
                if len(ssign2[(i,j)][1].keys()) == 0:
                    del ssign2[(i,j)]

                v[j]=v[j]+choice[inter][i][0]#the contributions from the directions of the local matching

                i=choice[inter][i][1]#moves according to local matching




        return v


###############################################################################



############# Some other (useful) functions. ##################################

def print_rank(answer):
    """Pretty printing of ranks of a bigraded homology theory.
    Input: a list of tuples((i,j),r) where (i,j) is the bigrading and r is the rank.
    """
    (t_min,t_max) = (min([x[0][0] for x in answer]),max([x[0][0] for x in answer]))
    (q_min,q_max) = (min([x[0][1] for x in answer]),max([x[0][1] for x in answer]))
    
    M = matrix(ZZ,(q_max-q_min)/2+2, t_max-t_min+2,sparse=True)
    
    for i in range(1,t_max-t_min+2):
        M[0,i]=i+t_min-1
    for i in range(1,(q_max-q_min)/2+2):
        M[i,0]=-2*i+q_max+2

    for ((t,q),r) in answer:
        M[(q_max-q)/2+1,t-t_min+1]=r

    #this is our output matrix. Now we will replace 0's by . (but not in top row or left column)
    splitting = M.str().split('\n')
    temp = splitting[0]
    ans = temp.replace('0','.',1)
    for i in range(1,(q_max-q_min)/2+2):
        temp = splitting[i]
        ind=temp.find(' 0')
        if ind != -1:
            almost = temp.split(' 0',1)
            stripped = almost[0].lstrip('[').lstrip()
            if stripped == '':
                temp=almost[0]+' 0'+almost[1].replace(' 0',' .')
            else:
                temp=almost[0]+' .'+almost[1].replace(' 0',' .')
        ans=ans+'\n'+temp
    print ans

def define_knot(name, print_progress = True):
    globals()[name]=KnotKh(eval(name+'_resolved_knot'), eval(name+'_cross'), print_progress)

def compute_knot(name, print_progress = False):
    globals()['Central_Knot']=KnotKh(eval(name+'_resolved_knot'), eval(name+'_cross'),print_progress)
    
"""
DOES NOT WORK WITH UNKNOT OR HOPF LINK COMPONENTS
"""
