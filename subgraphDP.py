## Subgraph Differential Privacy
## Implemented based on the pseudo code in:
## B. P. Nguyen, H. Ngo, J. Kim, and J. Kim, ‘‘Publishing graph data with subgraph differential privacy,’’ in Proc. Int. Workshop Inf. Secur. Appl. Jeju-do, South Korea: Springer, 2015, pp. 134–145.

## perturbGraph(G:nx.Graph, k = 5, epsilon = 5, beta=0.5):
## graph G, subgraph depth k, priv. budget epsilon, ratio of false to true edges beta

import numpy as np
import networkx as nx
import random
from scipy.stats import laplace

def selectEdges(G:nx.graph, k:int, beta):
    '''Select edges such that every k-vertex connected subgraph in G is perturbed.
    Done by doing Breadth First Search starting at v for v in G and selecting real and virtual edges.
    beta is the ratio of false edges to real edges, i.e. beta = |Ev| / |Er| -> |Ev| = |Er|*beta
    '''
    Es = []
    for v in G.nodes():
        Er_v = list(G.neighbors(v))
        Ev_v = [a for a in nx.bfs_tree(G, source=v, depth_limit=k).nodes() if a not in Er_v and a != v]
        sampleSize = min(int(len(Er_v)* beta), len(Ev_v)) ## prevent error if no false edges possible
        Es += [(min(v,node),max(v,node)) for node in Er_v + random.sample(Ev_v, sampleSize)]
    return list(set(Es))

def computeMetric(G:nx.Graph, Es):
    '''computes the normalized mutual friends metric of each pair of vertices in Es 
    m(i,j) = #mutualFriends / 2 * (1/d[i] + 1/d[j])
    Returns {(i,j): metric} for all edges in Es'''
    return {(i,j): len(set(G.neighbors(i)).intersection(G.neighbors(j)))/2 * (1/G.degree(i) + 1/G.degree(j)) for (i,j) in Es}

def scale(m:dict, alpha):
    '''changes the range of values in m from [0,1] to [alpha,1]'''
    ## multiplying by alpha gives range [0,alpha] so we want [0,1] -> [0,1-alpha] -> [alpha, 1]
    return {k: v*(1-alpha)+alpha for (k,v) in m.items()}

def addNoise(G:nx.Graph, Es, w, sigma, theta):
    '''
    graph G, edges Es, weights w, priv budg sigma, rewiring threshold theta
    random noise is injected into the weight of each e in Es
    for each e if w(e) < theta edge is rewired. if w(e) > theta stays the same
    random noise follows laplace with mu = 0 and sigma dep on privacy budget epsilon
    '''
    ## noise addition
    noise = laplace.rvs(loc = 0, scale = sigma, size = len(Es))
    newWeights = {Es[i]: w[Es[i]] + noise[i] for i in range(len(Es))}     ## edge : weight + noise
    ## graph rewiring
    Ganon = G.copy()
    for edge in Es:
        if newWeights[edge] <= theta:
            if edge in G.edges(): ## rewiring means edge changes state from real to virtual
                Ganon.remove_edge(*edge)
            else:
                Ganon.add_edge(*edge)
    return Ganon

def perturbGraph(G:nx.Graph, k = 5, epsilon = 5, beta=0.5):
    '''graph G, subgraph depth k, priv. budget epsilon, ratio of false to true edges beta'''
    Es = selectEdges(G,k, beta)   ## Es = set of edges of form (min, max) -> no repetitions
    Nk = k*(k-1)/2
    epsilon_i = epsilon/Nk
    sigma = -(1/np.log(2/(np.exp(epsilon_i)+1)))
    M = computeMetric(G, Es)    ## computes only for Es
    alpha = np.exp(-1/sigma)
    M = scale(M,alpha)
    w = {key: -sigma*np.log(value) for (key,value) in M.items()}        ## assigns a weight to each edge -> after metric scaling
    Gprime = addNoise(G, Es, w, sigma, theta = 0)
    return Gprime