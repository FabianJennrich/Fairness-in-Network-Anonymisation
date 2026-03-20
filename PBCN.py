## Privacy Preserving Approach Based on Clustering and Noise
## Implemented based on pseudo code in:
## H. Huang, D. Zhang, F. Xiao, K. Wang, J. Gu, and R. Wang, ‘‘Privacy-preserving approach PBCN in social network with differential privacy,’’ IEEE Trans. Netw. Service Manage., vol. 17, no. 2, pp. 931–945, Jun. 2020.

## PBCNanonymisation(G:nx.Graph)
## 10 groups,  epsilon = 0.1, n/100 as num changes

## NOTE: degree sequence for nx.havel_hakimi_graph is checked only for evenness, so sometimes fails due to non-graphical degree sequence when trying to anonymize very small graphs (100 nodes or less).

import numpy as np
import networkx as nx
import random
from sklearn.cluster import KMeans
from scipy.stats import dlaplace, binom

def groupConstruction(degSeq, T):
    '''
    Algorithm 1: ... based on K-means clustering
    Returns: dictionaries of clustering centers  C= clusterID: clusterCenter, groupNodes = nodeID: groupID, D = groupID:[nodes]
    '''
    ## dist(u,v) = |d(u) - d(v)|
    ## we assume degSeq is degree dictionary
    degSeqNodes = list(degSeq.keys())
    degSeqArr = np.array(list(degSeq.values())).reshape(-1,1)  ## create an array of degs while maintaining ordering
    mymeans = KMeans(n_clusters=T).fit(degSeqArr)

    unsortedC = {i: mymeans.cluster_centers_[i] for i in range(T)}  ## clusterID : cluster mean         clusterIDs are just 0-num clusters
    sortedClusterIDs = sorted(unsortedC, key= unsortedC.get, reverse=True)  ## clusterIDs by descending cluster mean
    ## oldID : newID 
    toNewID = {sortedClusterIDs[newid]: newid for newid in range(len(sortedClusterIDs))}

    C = {toNewID[oldID]: mymeans.cluster_centers_[oldID] for oldID in range(T)}
    nodeGroups = {degSeqNodes[i] : toNewID[mymeans.labels_[i]] for i in range(len(degSeqNodes))}          ## node: Old->New (old clusterID)
    D = {}
    for node in nodeGroups:
        if nodeGroups[node] in D:
            D[nodeGroups[node]].append(node)
        else:
            D[nodeGroups[node]] = [node]
    return C, nodeGroups, D

def randomDisturbance (G:nx.Graph, K, epsilon_1, nodeGroup, D):
    '''
    Algorithm 2: ... based on groups for pre-processing
    Returns reconstruction graph G'
    K is number of changed edges, nodeGroup is node: groupID, D is groupID : [nodes]
    '''
    T = len(D)
    Gprime = G.copy()
    epsilon_once = epsilon_1 / K
    for i in range(K):      ## do K times                                   ## 2
        (x,y) = random.choice(list(Gprime.edges()))                         ## a
        ## nodeGroup[x] nodeGroup[y] is i and j respectively
        ## Laplace noise eta_1 and eta_2 with 1/epsilon_once
        eta_1, eta_2 = dlaplace.rvs(1/epsilon_once, loc = 0, size = 2)   ## b
        iprime = (nodeGroup[x] + eta_1)%T
        jprime = (nodeGroup[y] + eta_2)%T       ## we use mod T so that we dont get overspill
        xprime = random.choice(D[iprime])       ## choose random nodes from new groups
        yprime = random.choice(D[jprime])
        if (xprime, yprime) not in Gprime.edges():
            Gprime.remove_edge(x, y)
            Gprime.add_edge(xprime,yprime)

    return Gprime

def noiseAllocation(Gprime:nx.Graph, D, C, epsilon_2):
    '''
    Algorithm 3: ... based on Laplace mechanism
    Returns noisy degree sequecne node: degree
    '''
    d = dict(Gprime.degree())
    m = Gprime.number_of_edges()

    Dprime = {i:[] for i in D}   ## same length as D, and all els = 0

    ## sensitivity S(f)_L(i) of f for ith group = Sf * sum Di / (2m) 
    ## C[i]*len(D[i]) = sum of Di degrees approx since we did perturbation
    Sf = {i: 2 * C[i]* len(D[i]) / m for i in C}       

    ## laplace noise eta_3 from Lap(S(f)_L(i) / epsilon_2)
    for i in C:     ## iterate through groupIDs
        Dinoise = dlaplace.rvs(Sf[i]/epsilon_2, loc=0, size = len(D[i]))
        ## for each node, we want to add noise
        Dprime[i] = [d[D[i][j]] + Dinoise[j] for j in range(len(D[i]))]
    dnoised = {D[i][j]: Dprime[i][j] for i in D for j in range(len(D[i]))}
    return dnoised

def buildNoiseEdge(Gprime:nx.Graph, noisedDegSeq):
    '''
    Algorithm 4: ... according to Havel Thm
    Returns reconstructed graph G'' (we call it Ganon for ease of reading)
    Takes G' (original graph with K changed edges), and noised degree sequence Dprime
    for (-) entries in Dprime, delete the edges constructed by Havel Theorem
    for (+) entries add edges constructed by Havel theorem
    NOTE: does not check if the degree sequences are graphical, makes them even by just adding 1 to deg of vertex 0, 
    if otherwise ungraphical just throws an error
    '''
    Ganon = Gprime.copy()
    ## we need an even sum of degrees to construct the degree sequence
    ## connect node of highest degree to other nodes of highest degree
    ## actually split into (+) and (-) parts, pass to havel_graph
    ## then add / remove edges
    Da = []
    Db = []
    Damap = {}
    Dbmap = {}
    for i in sorted(Gprime.nodes()):    ## we trust that these are in order 0-n     ## if they are not, there will be issues in Havel-Hakimi
        if noisedDegSeq[i] >= 0:
            Damap[len(Da)] = i
            Da.append(noisedDegSeq[i])
        else:   ## if < 0 i.e. (-)
            Dbmap[len(Db)] = i
            Db.append(-1*noisedDegSeq[i])
    if sum(Da) %2 ==1:
        Da[0] += 1
    if sum(Db) %2 ==1:
        Db[0] += 1
    DaHH = nx.havel_hakimi_graph(Da)
    DaHH = nx.relabel_nodes(DaHH, Damap, copy=False)
    DbHH = nx.havel_hakimi_graph(Db)
    DbHH = nx.relabel_nodes(DbHH, Dbmap, copy=False)        ## constructed edges based on havel thm, with intact node names

    Ganon.add_edges_from(DaHH.edges())
    for edge in DbHH.edges():
        if edge in Ganon.edges():
            Ganon.remove_edge(*edge)
    return Ganon

def postProcessing (Ganon:nx.Graph, epsilon_3):
    '''
    Algorithm 5: ... by adding node noises
    Returns reconstructed graph G'
    '''
    GanonPP = Ganon.copy()
    n = Ganon.number_of_nodes()
    d = [val for (_, val) in Ganon.degree()]
    d_noise = np.quantile(d, 0.05)  ## smallest 5% of nodes

    eta_4 = dlaplace.rvs(1/epsilon_3, loc=0, size = 1)[0]
    if eta_4 <0:
        smallDegs = [i for i in Ganon.nodes() if Ganon.degree(i) <= d_noise]
        if len(smallDegs) >= eta_4*-1:
            nodesToRemove = random.choices(smallDegs, k= eta_4*-1)
            GanonPP.remove_nodes_from(nodesToRemove)
        else:
            GanonPP.remove_nodes_from(smallDegs)
    if eta_4 >0:
        for i in range(eta_4):
            R = binom.rvs(n, p=d_noise/n, size = 1)[0] ## how many nodes to connect to
            nodesToConnect = random.choices(list(GanonPP.nodes()), k=R)
            edgesToAdd = [(n+i, node) for node in nodesToConnect]
            GanonPP.add_node(n+i)
            GanonPP.add_edges_from(edgesToAdd)
    return GanonPP

def PBCNanonymisation(G:nx.Graph):
    '''Produces anonymised graph
    T number of groups, K number of noise edges, epsilon privacy budget
    smaller epsilon the more edges are added -> longer run time
    Returns: Ganon anonymous nx.Graph
    '''
    epsilon=0.1
    K = int(G.number_of_edges()/100)+1
    T=10

    degSeq = dict(G.degree())
    C, groupNodes, D = groupConstruction(degSeq=degSeq, T=T)
    Gprime = randomDisturbance(G, K, epsilon_1=epsilon/3, nodeGroup=groupNodes, D=D)
    degSeqNoised = noiseAllocation(Gprime, D, C, epsilon/3)
    Ganon = buildNoiseEdge(Gprime, degSeqNoised)
    Ganon = postProcessing(Ganon, epsilon/3)    ## d_noise is smallest 5% of degrees
    return Ganon