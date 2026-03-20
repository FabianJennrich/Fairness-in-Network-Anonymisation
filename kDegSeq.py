## k-Degree Sequence Anonymisation
## Implemented based on the description in:
## J. Casas-Roma, J. Herrera-Joancomartí, and V. Torra, ‘‘K-degree anonymity and edge selection: Improving data utility in large networks,’’ Knowl. Inf. Syst., vol. 50, no. 2, pp. 447–474, Feb. 2017.

## k_degSeq_Anon(G:nx.Graph(), k= 10)
## creates k=10 degree anonymity

import networkx as nx
import numpy as np
from microagg1d import univariate_microaggregation
import random

def greedyDegSeqAnon(m, probs, maxIt = 20):
    ''' Performs greedy optimization for target degrees (mj). Values for mj are selected according to prob distribution based on size of mj1 and mj2. The process is finished when a solution is found s.t. sum_{j=1 to p} of (mj) == 0, 
    or when we go a certain number (20) of iterations without change.'''
    ## just get some option to start
    sumM = 1 ## do the definition at least once
    while sumM % 2 == 1:
        bestM = {j: m[j][0] if random.random() <= probs[j][0] else m[j][1] for j in m}
        sumM = abs(sum(bestM.values()))

    ## select new option
    daysSinceLastChange = 0 ## iteration counter
    while daysSinceLastChange < maxIt:
        newM = {j: m[j][0] if random.random() <= probs[j][0] else m[j][1] for j in m}
        if abs(sum(newM.values())) % 2 == 0 and abs(sum(newM.values())) < sumM:
            bestM = newM
            sumM = abs(sum(newM.values()))
            daysSinceLastChange = 0
        else:
            daysSinceLastChange += 1
    return bestM

def k_degSeq_Anon(G:nx.Graph(), k= 10):
    ''' Anonymizes the graph G by k-degree sequence anonymization. The degree sequence is partitioned using UMGA. Each group is assigned a target degree using greedy optimization of all degree changes. Edges are added, removed, and switched to achieve the target degrees.
    Returns: Ganon anonymous nx.Graph'''
    ## clustering the degree sequence
    d = sorted(dict(G.degree()).items(), key = lambda x: x[1], reverse = True)
    nodesL = [n for (n,_) in d]
    degL = [deg for (_,deg) in d]
    clusters = univariate_microaggregation(degL, k=k)
    ## list of cluster numbers, ordered as in nodesL

    ## want to create mj1 and mj2 for each group j
    ## need the avg deg of each cluster
    ## create dictionary s.t. key is cluster number, and value is list of corresponding nodes
    clusters_degsDict = {j:[] for j in set(clusters)} ## get only the cluster numbers
    for i in range(len(degL)):
        clusters_degsDict[clusters[i]].append(degL[i])
    ## now we have clustersDict as {cluster number: [di,...,dn]} for i,...,n in cluster
    sumd = {j:(sum(clusters_degsDict[j]), len(clusters_degsDict[j])) for j in clusters_degsDict}
    m = {j:(sumd[j][0]-sumd[j][1]*np.floor(sumd[j][0] / sumd[j][1]), sumd[j][0]-sumd[j][1]*np.ceil(sumd[j][0] / sumd[j][1])) for j in sumd}
    probs = {j:(1- m[j][0]/(m[j][0]-m[j][1]),1 + m[j][1]/(m[j][0]-m[j][1])) for j in m}
    ## use implied def for probs rather than stated (stated def not in range[0,1])
    ## sometimes encounters div 0 error, but that seems fine, since => mj1=mj2=0, so no change necc
    
    bestM = greedyDegSeqAnon(m, probs)
    ## bestM is a dictionary {j: mj1 or mj2}
    ## we know which one based on the sign, if (-) then mj2
    targets = {j: np.floor(sumd[j][0]/sumd[j][1]) if bestM[j]>0 else np.ceil(sumd[j][0]/sumd[j][1]) for j in bestM}
    
    delta = {nodesL[i]:targets[clusters[i]]-degL[i] for i in range(len(clusters))}
    delta_min = [node for (node,defc) in delta.items() if defc<0]
    delta_plus = [node for (node,defc) in delta.items() if defc>0]

    Ganon = G.copy()
    failureCount = 0 ## ex: delta+ = [vi,vj], and (vi,vj) in E -> no edge addition is possible
    ## if too many failures occur, returns graph to date
    
    ## we do random edge selection 
    while sum(delta.values()) < 0 and failureCount < 50: ## edge deletion
        try:
            vi,vj = random.sample(delta_min, k=2)
            vk = random.choice([vertex for vertex in Ganon.neighbors(vi) if vertex != vj])
            vl = random.choice([vertex for vertex in Ganon.neighbors(vj) if vertex != vi and vertex not in Ganon.neighbors(vk)])
            Ganon.remove_edges_from([(vi,vk),(vj,vl)])
            Ganon.add_edge(vk,vl)
            delta[vi] += 1 ## - -1
            delta[vj] += 1
            delta_min = [node for (node,defc) in delta.items() if defc<0]
            failureCount=0
        except:
            failureCount += 1
        
    while sum(delta.values()) > 0 and failureCount < 20: ## edge addition
        try:
            vi,vj = random.sample(delta_plus, k=2)
            while (vi,vj) in Ganon.edges():
                vi,vj = random.sample(delta_plus, k=2)
            Ganon.add_edge(vi,vj)
            delta[vi] -= 1
            delta[vj] -= 1
            delta_plus = [node for (node,defc) in delta.items() if defc>0]    ## edge switch
            failureCount = 0
        except:
            failureCount += 1
    ## while there are nodes with not the target degree:
    while sum(np.abs(list(delta.values()))) != 0 and failureCount <20: 
        try:
            vi = random.choice(delta_min)
            vj = random.choice(delta_plus)
            vk = random.choice([v for v in Ganon.neighbors(vi) if v not in Ganon.neighbors(vj) and v != vj])
            Ganon.remove_edge(vi,vk)
            Ganon.add_edge(vj,vk)
            delta[vi] += 1
            delta[vj] -= 1
            delta_min = [node for (node,defc) in delta.items() if defc<0]
            delta_plus = [node for (node,defc) in delta.items() if defc>0]
            failureCount = 0
        except:
            failureCount += 1
    return Ganon