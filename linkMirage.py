## Link Mirage Anonymisation
## Implemented based on pseudo code in:
## C. Liu and P. Mittal, “Linkmirage: Enabling privacy-preserving analytics on social relationships,” in Proc. NDSS, San Diego, CA, USA, 2016, pp. 1–15.

## linkMirage(G:nx.Graph, t=5)
## t the random walk length

import networkx as nx
import random
from itertools import combinations

def randomWalk(G:nx.Graph, u, distance):
    '''Returns the terminal vertex of a random walk of length distance starting at u'''
    vertex = u
    for i in range(distance):
        vertex = random.choice(list(G.neighbors(vertex)))
    return vertex

def transform(G:nx.Graph, t=5, M=10):
    '''Returns a list of perturbed edges for graph G
    t is the random walk length -> chose 5 from paper results
    M is the number of tries to find a suitable edge'''
    GprimeEdges = []
    for u in G.nodes(): ## foreach u in G
        count = 1
        for v in G.neighbors(u):    ## foreach neighbor v of u
            loop = 1
            while True:
                z = randomWalk(G, v, t-1)   ## t-1 hop rand walk from v
                loop += 1
                if loop <= M:
                    break
                if u != z and (u,z) not in GprimeEdges and (z,u) not in GprimeEdges:
                    break
            if loop <= M:
                if G.degree(u) > 1:
                    prob = (0.5 * G.degree(u) - 1) / (G.degree(u) - 1)
                else:
                    prob = 0
                if count == 1:
                    GprimeEdges.append((u,z))
                elif random.random() <= prob:
                    GprimeEdges.append((u,z))
                count += 1      ## increase count only if an edge was added -> elsewise move on to the next neighbor without increasing count
    return GprimeEdges

def staticPerturbation(G:nx.Graph, Gprime:nx.Graph, C, t):
    '''Perturbs the changed communities as in linkPerturbation 
    takes original graph G, graph to be modified Gprime, and list of node sets C
    modifies in place Gprime'''
    ## do transform for subgraphs induced by the communities
    for nodes in C:
        subGraph = G.subgraph(nodes)
        newEdges = transform(subGraph, t=t)
        if newEdges:    ## if newEdges not empty
            Gprime.add_edges_from(newEdges) ## vertices should stay the same because we do everything by node id

def marginalNodes(G:nx.Graph, a, b):
    '''Returns lists of marginal nodes in communities a and b, where a and b are sets of nodes in graph G
    va = [(node, num neighbors in b)], vb = [(node, num neighbors in a)]'''
    ## nodes in a s.t. the intersection of its neighbors with b is nonempty
    return [(v, len(set(G.neighbors(v)).intersection(b))) for v in a if set(G.neighbors(v)).intersection(b)], [(v,len(set(G.neighbors(v)).intersection(a))) for v in b if set(G.neighbors(v)).intersection(a)]

def linkMirage(G:nx.Graph, t=5):
    '''Anonymizes G using LinkMirage algorithm with t the random walk length used.
    Returns: Gprime anonymous nx.Graph'''
    Gprime = nx.Graph()
    Gprime.add_nodes_from(G.nodes())    ## has same nodes
    ## if t == 0:
    C0 = nx.community.greedy_modularity_communities(G)  ## cluster G0 to get C0     ## C0 is a list of frozen sets of nodes -> immutable set
    C0_ch = C0  ## label C0 as changed, i.e. C0_ch = C0

    ## Dynamic Clustering -> unnecc

    ## Selective Perturbation
    ## get changed and unchanged communities -> all changed so just C0_ch
    ## G't_un = G'(t-1)_un -> transfer unchanged communities as is

    ## intra-cluster Perturbation
    staticPerturbation(G, Gprime, C0_ch, t)    ## perturb Ct_ch for G't_ch by static method    ## needs G for perturbation, modify Gprime in place

    ## inter-cluster perturbation
    for (a,b) in combinations(C0_ch, 2):    ## all combos of communities
        ## if both unchanged -> unnecc, all changed
        va, vb = marginalNodes(G,a,b)   ## get marginal node sets
        len_E_ab = sum(deg for _,deg in va) ## total degree of marginal nodes va
        for (i,di) in va:
            for (j,dj) in vb:
                prob = di*dj*len(va) / (len_E_ab * (len(va)+len(vb)))
                if random.random() <= prob:
                    Gprime.add_edge(i,j)
    ## end Selective Perturbation
    return Gprime
