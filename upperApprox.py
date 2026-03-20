## Upper Approximation based Anonymisation
## Implemented based on working example in:
## S. Kumar and P. Kumar, ’’Upper approximation based privacy preserving in online social networks,’’ Expert Syst. Appl., vol. 88, pp. 276–289, Dec. 2017.

## PPGPworkingEx(G:nx.Graph, delta = 0.2)
## delta is the threshold for minimum Jaccard Similarity between nodes

import numpy as np
import networkx as nx
from itertools import combinations

def d_neighborhood(G, d):
    '''
    Returns a dictionary with as keys the nodes, and as values the d-neighborhood of that node as a set
    time complexity of this should be O(n*h^d) where h avg node degree, and d neighborhood size
    takes about 1 sec for d=3,4 h = 50, n = 1000, about 10 secs for n=10k d=2, h=100
    '''
    DDS = {}
    for v in G.nodes(): ## gets the d neighborhood of all nodes in G saves in a dict
        n0 = {v}
        totalset_v = set()
        for i in range(d):
            n1 = set()
            for node in n0:
                n1.update(G.neighbors(node))
            totalset_v.update(n0)
            n0 = n1 - totalset_v
        totalset_v.update(n0)
        DDS[v] = totalset_v
    return DDS

def PPGPworkingEx(G:nx.Graph, delta = 0.2):
    '''
    Returns the anonymized graph based on the implied definitions from the working example
    that is the KiteGraph example
    '''
    Sbar = d_neighborhood(G, 2)
    JS = {vi:{vj: len(set(G.neighbors(vi)).intersection(G.neighbors(vj))) / len(set(G.neighbors(vi)).union(G.neighbors(vj))) for vj in G.nodes()} for vi in G.nodes()}
    Shat = {}
    for vi in G.nodes():
        myset_vi = set()
        for vl in Sbar[vi]:
            myset_vi.update([vj for vj in Sbar[vl] if JS[vi][vj] >= delta])
        Shat[vi] = myset_vi
    anonG = nx.Graph()
    anonG.add_nodes_from(G.nodes())
    for vi in G.nodes():
        for vj in Shat[vi]:
            if vi != vj:
                anonG.add_edge(vi,vj)
    return anonG