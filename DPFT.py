## Differential Privacy with Field Theory
## Implemented based on pseudo code in:
## DPFT (Differential Privacy with Field Theory (DPFT)): H. Zhu, X. Zuo, and M. Xie, ‘‘DP-FT: A differential privacy graph gen- eration with field theory for social network data release,’’ IEEE Access, vol. 7, pp. 164304–164319, 2019.

## DP_FT(G:nx.Graph, epsilon = 1, t:int = 5):
## epsilon is noise level, t is random walk length
## Line number comments refer to lines in the pseudo code of the DPFT paper

import numpy as np
import networkx as nx
import itertools
import random
from scipy.stats import dlaplace

def distance_compute(G:nx.Graph, t:int):
    ## 1: Initialize TransMatrix and SimMatrix where TransMatrix has size |V|x|V| and 
    d = G.degree()
    SimDict = {i:{j:1/d[i] if (i,j) in G.edges() else 0 for j in G.nodes()} for i in G.nodes()}
    SimMatrix = pd.DataFrame(SimDict)   ## for the sake of my sanity we use pandas (ease of use)
    TemMatrix = SimMatrix.copy()
    HelperMatrix = TemMatrix.copy()
    TransMatrix = SimMatrix.copy()      ## does not change from here
    walksteps = 1
    while walksteps < t:        ## 2: repeat 5: until walksteps >= t
        ## 3: random walk calculate 2-step transistion prob vi to vi
        ## 2-step transition probability is given by prob that the m+2th stop is j given that the mth stop is i
        ## prob that i goes to k * prob k goes to j; sum of s[i][k]*s[k][j] for all  k
        for i in G.nodes():
            for j in G.nodes():     ## using pseudocode method, there is interference with the updating columns
                ## df["col"][row_indexer] = value Use `df.loc[row_indexer, "col"] = values` instead,
                ##TemMatrix[i][j] = sum(HelperMatrix[i][:]*TransMatrix[:][j])
                ## sum is over intersection of neighbors of vi and vj so:
                vk = [ni for ni in G.neighbors(i) if ni in G.neighbors(j)]
                TemMatrix.loc[j,i] = sum((HelperMatrix[i][:]*TransMatrix.T[j][:])[vk])
        ## 4: Update SimMatrix
        SimMatrix = SimMatrix + TemMatrix
        HelperMatrix = TemMatrix.copy()
        walksteps += 1
    ## 6:
    return SimMatrix

def fast_distance_compute(G:nx.Graph, t:int):
    '''Compute the SimMatrix of the graph G for a random walk length of t
    '''
    ## sim[i][j] = pij + pij^2 + pij ^3 + ... + pij^t where pij^n is the probability to go from i to j in exactly n steps
    ## we can compute this by computing powers of the transition matrix, (to get pij^n) and adding them
    d = dict(G.degree())
    T = [[(1/d[i])*((i,j) in G.edges()) for i in G.nodes()] for j in G.nodes()]
    T = np.array(T)
    ## getting the transition matrix. We assume in doing this that the nodes are in a sensible order
    helper = T.copy()
    SimMatrix = np.zeros((G.number_of_nodes(),G.number_of_nodes()))
    for i in range(1,t+1):
        ## helper = T^i
        SimMatrix = np.add(SimMatrix,  helper) ## SimMatrix = sum of those
        helper = np.matmul(helper,T)  ## helper = T^i+1
    return SimMatrix

def DP_FT(G:nx.Graph, epsilon = 1, t:int = 5):
    '''
    Input: graph G, privacy budget epsilon, number of steps for random walk t (default 2 as in original paper)
    Output: list of edges E_T
    note: original paper uses epsilon in [0.5, 3] in steps of 0.5, they got good results for epsilon <= 1
    '''
    ## 1: initialization of SimMatrix, M_force, E_T, and queue
    
    F = {i:{j:0 for j in G.nodes()} for i in G.nodes()}
    E_T = []
    queue = [] 
    V = list(G.nodes())     ## 2      ## list form so that we can sensibly modify it later
    d = dict(G.degree())    ## 3    ## we call it d rather than S so that getting the degree of a node is more intuitive
    SimMatrix = fast_distance_compute(G,t)   ## 4, 1 (initialization)
    ## 5,6,7: for each node pair vi vj in G, NodeSim(i,j) = SimMatrix[i][j] -> just use SimMatrix[i,j]
    candidates = {}         ## vi_candidates is used again later so we initialize it outside the loop

    ## ensures things still work correctly when nodes are arbitrary names not 0-n
    ## node Translation Dictionary
    nTD = dict(zip(sorted(G.nodes()), range(G.number_of_nodes())))
    
    for vi in V:            ## 8, 16
        candidates[vi] = [vv for vv in V if SimMatrix[nTD[vi]][nTD[vv]] != 0 and vi != vv] ## 9   ## seems like we prob want vi = vv excluded
        ## 10,14: if vi_cand not empty then:
        for vv in candidates[vi]:    ## 11, 13
            ##12: M_force[i][v] = F(i,v)    where F(a,b) = d(a) * d(b) * NodeSim(a,b)^n where NodeSim(i,j) = SimMatrix[i][j]
            ##       we omit ^n because we use n = 1, as in experiments from orig. paper
            F[vi][vv] = d[vi]*d[vv]*SimMatrix[nTD[vi]][nTD[vv]]   ## vertex indexing is the same as the vertices by construction
    
    ## 16: // adding Laplacian noise to degree sequence S to get S_tilde -> modified in place, S is called d, for ease of use
    if epsilon != 0:        ## 17, 22
        for vi in d.keys(): ## 18, 21   ## d is node:degree dictionary
            di_tilde = min(len(V)-1, max(1, d[vi] + dlaplace.rvs(2/epsilon)))    ## 19
            ## should probably be len(V)-1 since we cannot connect to self
            d[vi] = di_tilde## 20
    
    while len(V) > 0:       ## 23, 45
        ## 24: // Mechanism for preferentially selecting a node
        ## 25: sample vi in V with prob di/sum(dk, k in V)
        prob = [d[vi] for vi in V]   ## 25
        vi = random.choices(population=V, weights=prob)[0]      ## 25
        if len(candidates[vi]) != 0: ## 26, 35      ## if vi has potential connections
            ## 27: // Mechanism for preferentially generating an edge
            ## 28: sample vj in vi_candidate with prob p(eij|vi) = F(i,j)/sum(F(i,k), for k in vi_candidates)
            sumFik = sum([F[vi][k] for k in candidates[vi]]) ## 28
            prob = [F[vi][vj]/sumFik for vj in candidates[vi]]  ## 28
            vj = random.choices(population=candidates[vi], weights=prob)[0]   ## 28
            candidates[vi].remove(vj)   ## 29
            candidates[vj].remove(vi)   ## 29
            if ((vi, vj) not in E_T) and ((vj,vi) not in E_T):  ## 30, 32   ## add the edge if it doesnt exist yet
                E_T.append((vi,vj)) ## 31
                d[vi] -= 1          ## 31
                d[vj] -= 1          ## 31
            if d[vj] == 0:   ## 39      ## remove vj from V if vj has no potential connections left (this is in here bc we do this only if vj has been selected this round)
                V.remove(vj)
        else:   ## 33, 35               ## if there are no candidates for vi, add to q
            queue.append(vi)

        if d[vi] == 0:   ## 36          ## remove vi from V if vi has no potential connections left
            V.remove(vi)
        if set(queue + V) == set(queue):    ## 41
            break   ## exit the while loop if all vertices in V are already in the queue
    ## 45: end while

    while len(queue) > 0:       ## 46
        vi = random.choice(queue)         ## 47   ## selection method is uniform at random
        queue.remove(vi)    ## we have now addressed this elm. and all other bits are done from V (not mentioned in pseudocode, implied by naming of queue)
        ## generate list of unused pot. edges -> check for all potential edges if they have been used
        potential_vj = [v for v in V if v!=vi]
        for vj in potential_vj:     
            if (vi,vj) in E_T or (vj,vi) in E_T:
                potential_vj.remove(vj)

        if len(potential_vj) > 0:            ## 48: if V not empty   ## pseudocode would run forever if all pot. edges used up
            workingV = V
            sumDegk = np.sum([d[k] for k in V])
            prob = [d[vj]/sumDegk for vj in V]
            vj = random.choices(population=V, weights= prob)[0]     ## 50
            if vj == vi or (vi,vj) in E_T or (vj,vi) in E_T:     ## 49,51    ## wait until none of these are true
                workingV.remove(vj)
                sumDegk = np.sum([d[k] for k in workingV])
                prob = [d[vj]/sumDegk for vj in workingV]
                vj = random.choices(population=workingV, weights= prob)[0] ## 50
            d[vj] -= 1          ## 52
            E_T.append((vi,vj)) ## 52
            if d[vj] == 0:   ## 53, 55
                V.remove(vj)    ## 54
        else:   ## if V empty / if there are no possible edges for vi
            if vi in V: ## else we dont have to bc its already gone
                V.remove(vi)    ## remove it from nodeslist
    
    ## 58: return E_T
    anonG = nx.Graph()
    anonG.add_nodes_from(G.nodes())
    anonG.add_edges_from(E_T)
    return anonG    ## this algorithm returns an anonymized graph rather than an edgelist