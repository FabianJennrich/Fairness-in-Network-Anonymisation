## K- Degree anonymisation with Vertex and Edge Modification
## Implemented based on pseudo code in:
## T.Ma,Y.Zhang,J.Cao,J.Shen,M.Tang,Y.Tian,A.Al-Dhelaan,and M. Al-Rodhaan, ‘‘KDVEM: A k-degree anonymity with vertex and edge modification algorithm,’’ Computing, vol. 97, no. 12, pp. 1165–1184, 2015.

## KDVEMapproach(G:nx.Graph, k = 10)
## creates k=10 degree anonymity

import networkx as nx
import random

def greedy_partition(degDict, k):
    '''Inputs: degDict = {node:degree}, anon level k, 
    Rerturns: targetDegDict = node: targetDeg, and dprime = targetDeg : [nodes]'''
    d = sorted(degDict.items(), key = lambda item : item[1], reverse=True)
    ## dprime = {target degree: [nodes]}, initialized with maxDeg : [k largest deg nodes]
    dprime = {d[0][1]:[pair[0] for pair in d[0:k]]}    
    target_d = d[0][1]
    i = k
    ## iterate through the degree pairs, except the last k (they are necc. added to a group)
    while i < len(d)-k:  ## iterate through the degree pairs, except the last k (they are necc. added to a group)
        (vi, dvi) = d[i]
        ## compute costs to merge and to make new group
        ## since we do addition to make series same, cost_merge= maxDeg_i-deg v and cost_new = k*deg_v - sum(deg_v+1,...,deg_v+k)
        C_merge = target_d - dvi
        C_new = k*dvi - sum([pair[1] for pair in d[i:i+k]]) ## sum of next k degrees
        if C_merge <= C_new:        ## default has to be merge otherwise there are issues
            dprime[target_d].append(vi)
            i+= 1
        else:
            target_d = dvi
            ## get next k nodes in deg sequence and assign to dprime with key= current deg
            dprime[target_d] = [pair[0] for pair in d[i:i+k]]    
            i += k 
    dprime[target_d] = dprime[target_d] + [pair[0] for pair in d[i:]]
    targetDegDict = {node: targetDeg for targetDeg in dprime for node in dprime[targetDeg]}
    return targetDegDict, dprime

def find_candidates(G:nx.Graph, Vplus, id_communities):
    '''     Creates a list of candidates for each vertex, by iterating through vertex community by increasing level (with full graph at top level), and selecting edges in order of increasing social distance between nodes.
    INPUT: G graph, degDef = {node: degree deficiency}
    Returns: candidatesDict node: [(candidate, dist cand to node)] which contains only nodes with degree deficiency
    '''
    n = G.number_of_nodes()
    candidates_Dict = {}
    for v in Vplus:
        candidates = []
        communities = id_communities[v]
        for lvl in range(len(communities)):
            temp = []
            for u in Vplus:
                if u!=v and id_communities[u][lvl] == communities[lvl] and (u,v) not in G.edges():
                    try:
                        distance = nx.shortest_path_length(G,u,v)
                    except:
                        distance = n
                    temp.append((u,distance))
                ## end if
            ## end for
            temp = sorted(temp, key = lambda item:item[1])
            candidates = candidates + temp
        ## end for
        candidates_Dict[v] = candidates
    return candidates_Dict      

def community_multilevel(G:nx.Graph):
    ''' Creates the community partition for the nodes at several levels using louvain partition.
    Returns: id_communities = {nodeID:[cID_l0,cID_l1, ... ,cID_lmax]}
    dictionary with nodes as keys, value is a list that are the commIDs at different levels
    levels sorted lowest to highest
    '''
    partitions = nx.community.louvain_partitions(G)
    id_comm = {nID:[] for nID in G.nodes()}
    for level in partitions: ## lowest to highest level
        lvlDict = {node:pID for pID in range(len(level)) for node in level[pID]}
        id_comm_helper = {nID: id_comm[nID] + [lvlDict[nID]] for nID in G.nodes()}
        id_comm = id_comm_helper
    ## add graph as highest lvl comm
    id_communities = {nID: id_comm[nID]+ [0] for nID in G.nodes()} 
    return id_communities
    
def graph_modification(G:nx.Graph, k, d, d_target, dprime, id_communities):
    '''
    Transforms G to its k-anonymous graph Gprime. For each v in Vplus, first try to increase its degree by connecting it with candidates that exists in the original graph. Candidates c fulfil that c needs degree increase, and (c,v) not in original graph. After modifying the graph through adding edges, if there exists vertices still need to increase its degree, we try to modify the graph through adding vertices to the graph as in [pseudeoNode anon]
    Returns: Gprime nx.Graph
    '''    
    Gprime = G.copy()
    degSeq = sorted(d.items(), key = lambda item : item[1], reverse=True)  ## list of (node:degree)
    degDeficiency = {node: d_target[node]-d[node] for node in G.nodes()}
    Vplus = [node for node in degDeficiency if degDeficiency[node] >0]
    candidates = find_candidates(G, Vplus, id_communities)                  ## 11
    for v in Vplus:                                                 ## 9
        temp = degDeficiency[v]
        while temp > 0:                                             ## 13
            if len(candidates[v]) <= 0:                             ## 18
                break   ## exit while loop
            candidate = candidates[v].pop(0)    ## (node, distance) ## 14
            ## if candidate != null:
            if degDeficiency[candidate[0]] > 0 and (candidate[0],v) not in Gprime.edges():
                Gprime.add_edge(v, candidate[0])
                temp -= 1
                degDeficiency[candidate[0]] -= 1
        ## end while                                                ## 21
        if temp == 0:                                               ## 22
            Vplus = [node for node in Vplus if node != v]
        degDeficiency[v] = temp
    ## end for
    
    ## addVertex as in pseudoNode anonymization
    ## Chester S et al (2013) Why Waldo befriended the dummy? k-Anonymization of social networks with pseudo-nodes. Soc Netw Anal Min 3(3):381–399
    
    ## add m = (1+max(md,k))(mod2) + max(md,k)     md = di - dj where di largest deg, dj smallest deg, within any partition -> so largest val in deg Def
    md = max(degDeficiency.values())
    m = (1+max(md,k))%2 + max(md,k)
    newNodes = list(range(Gprime.number_of_nodes()+1,Gprime.number_of_nodes()+m+1))
    Gprime.add_nodes_from(newNodes)

    ## iterate through new nodes, and add edge to old nodes -> cycle through these slower
    newNodeID = 0
    for node in degDeficiency:
        for i in range(degDeficiency[node]):
            Gprime.add_edge(node, newNodes[newNodeID%len(newNodes)])
            newNodeID +=1
    d = int(newNodeID/len(newNodes)) + 1
    ## d is max degree of new nodes
    
    if (d in d_target.values() and d-1 in d_target.values()) or newNodeID%len(newNodes) == 0: ## also k-anon if all new nodes have same deg
        return Gprime
        
    ## vertices with deg d-1 pair them and add an edge between them
    leftoverNodes = newNodes[newNodeID%len(newNodes):]      ## list of length m-td mod m with nodes of degree d-1
    random.shuffle(leftoverNodes)
    for i in range(0,len(leftoverNodes)-1,2):
        Gprime.add_edge(leftoverNodes[i],leftoverNodes[i+1])
    ## that is [newNodeID to end of list] of new nodes -> if m-td mod m is even and m >= k then we are done
    if len(leftoverNodes)%2 ==0:
        return Gprime
    ## if m-td mod d is odd, then we leave out 1 vertex -> m-1 even and >=k 
    ## add an edge from r to two other new vertices (from full list) -> have degree
    edgesDeg_d = newNodes[0:newNodeID%len(newNodes)] + leftoverNodes[:-1]   ## all but the last element
    random.shuffle(edgesDeg_d)
    Gprime.add_edges_from([(leftoverNodes[-1],edgesDeg_d[0]),(leftoverNodes[-1],edgesDeg_d[1])])
    for i in range(2,len(edgesDeg_d)-1,2):  ## pair off, but els 0 and 1 are already taken
        Gprime.add_edge(edgesDeg_d[i],edgesDeg_d[i+1])
    return Gprime

def KDVEMapproach(G:nx.Graph, k = 10):
    '''
    Anonymizes G using KDVEM to create k-degree anonymous graph.
    '''
    d = dict(G.degree())
    id_communities = community_multilevel(G)
    target_deg_dict, dprime = greedy_partition(d, k)
    Ganon = graph_modification(G, k, d, target_deg_dict, dprime, id_communities)
    return Ganon