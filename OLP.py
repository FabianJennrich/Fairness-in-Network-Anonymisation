## Optimisation by Link Perturbation
## Implemented based on the pseudocode given in:
## T. Wu, G. Ming, X. Xian, W. Wang, S. Qiao, and G. Xu, ‘‘Structural predictability optimization against inference attacks in data publishing,’’ IEEE Access, vol. 7, pp. 92119–92136, 2019.

import networkx as nx
import scipy
import random
import numpy as np
from collections import Counter

def linkImportanceMeasuring(G:nx.Graph, c = 0.1, alpha = 0.2):
    '''
    Computes list of edges ranked by ascending importance
    Inputs: G graph to anonymize
        c float parameter for random walks
        alpha used to compute number of random walks
    Output: Rlist ranked edge list in ascending order of importance
    '''
    n = G.number_of_nodes()
    N = alpha * n
    maxIter = 100
    
    Pt = np.full(n, 1/n) ## P is the distribution to choose from, Pt = P0
    ## Pt+1 = (1-c)S^T Pt + c/n * 1 where 1 a suitable 1 vector, S = DA s.t. D=1/di 
    ## is the probability of ending a RW on a node at time t -> no starting node given
    A = nx.adjacency_matrix(G)
    D = scipy.sparse.diags_array([1/d for d in dict(G.degree()).values()])
    S = D@A
    covern_1 = np.full(n, c/n)  

    Rt = random.choices(range(n), weights = Pt, k = n)

    W = []
    for t in range(maxIter):
        Ptplus1 = (1-c)* S.T@Pt + covern_1
        Rtplus1 = random.choices(range(n), weights = Ptplus1, k = n)
        W += [i for i in zip(Rt, Rtplus1)]
        Pt = Ptplus1
        Rt = Rtplus1
    ## based on my understanding of the paper, we sample at each step t from the distribution
    ## we do not consider conditional probability

    W = [(min(u,v), max(u,v)) for (u,v) in W]
    ## update link importance matrix Q based on W, and get list Rlist in order
    c = Counter(W)
    Rlist = sorted(list(G.edges()), key = lambda e: c[e])
    return Rlist

def optimizationLinkPerturbation(G:nx.Graph, beta = 0.1):
    '''Takes G and perturbation ratio beta, and computes a graph with only the (1-beta) most important edges
    Inputs: G graph to anonymize
        beta fraction of edges removed
    Output: Gstar nx.Graph the anonymized graph
    '''
    Rlist = linkImportanceMeasuring(G)
    Gstar = G.copy()
    Gstar.remove_edges_from(Rlist[:int(beta*G.number_of_edges())])
    return Gstar
