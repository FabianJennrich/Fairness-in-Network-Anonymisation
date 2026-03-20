## Pseudo Node Anonymisation
## Implemented based on the pseudo code in:
## Chester S et al (2013) Why Waldo befriended the dummy? k-Anonymization of social networks with pseudo-nodes. Soc Netw Anal Min 3(3):381–399

## pseudoNodeAnon(G:nx.Graph, k=None, doGreedy=False)
## default for k is 10

## Note: for graphs with many nodes, computing the degree sequence recursively causes max recursion depth error, in which case greedy partitioning is used

import numpy as np
import networkx as nx
import random

def Start_unspec(degSeq,start,cost, k, x):
    ''' rightmost partition is given by [Start(x), x], next is given by [Start(Start(x)-1), Start(x)-1]
    where x is an index position
    '''
    if start[x] == None:
        if x < 2*k:
            start[x] =  1    ## bc we have constructed the indexing s.t. 0 should not be reached
        else:
            start[x] = Pos_Split_unspec(degSeq,start,cost, k, x)
    return start[x]

def Pos_Split_unspec(degSeq,start,cost, k, x):
    '''argmin over i in [max(k,x-2k+1), x-k]    of max(Cost(1,i-1), delta(i,x))
    assume that argmin returns maximal point at which a function is minimised (i.e. largest index)
    NOTE: condiditon to enter is x>=2k'''
    myiterator = list(range(x-k, max(k,x- 2*k +1)-1, -1))  ## we go in reverse order, so that we get maximal pt where func is min.
    ## want i, not the actual argmin so:
    return myiterator[np.argmin([max(Cost_unspecified(degSeq,start,cost, k, 1, i-1), delta_unspecified(degSeq, i, x)) for i in myiterator])]

def delta_unspecified(degSeq, x, y):
    '''return dx - dy
    where x and y are positions in the deg sequence, counting from the right
    degSeq = [(node,maxdeg),...,(node,mindeg)]'''
    return degSeq[x][1] - degSeq[y][1]  ## there will be an error if it ever tries to use the 0 index

def Cost_unspecified(degSeq,start,cost, k, _, x):
    '''delta(1,x) if x<2k, Cost_Split if x>=2k'''
    if cost[x] == None:
        if x < 2*k:
            cost[x] = delta_unspecified(degSeq, 1, x)
        else:
            cost[x] = Cost_Split_unspec(degSeq,start,cost, k, x)
    return cost[x]

def Cost_Split_unspec(degSeq,start,cost, k, x):
    '''minimum over i in [max(k, x-2k+1), x-k]      of (max(Cost(1,i-1), delta(i,x)))
    NOTE: condiditon for entering is k >= 2k'''
    myiterator = list(range(max(k, x-2*k+1), x-k+1)) ## want to include x-k so we have to use x-k+1
    return min([max(Cost_unspecified(degSeq,start,cost, k, 1, i-1), delta_unspecified(degSeq, i, x)) for i in myiterator])

def greedyGrouping(degSeq, k):
    '''a greedy method for determining the partition of deg seq
    this should allow it to run with more vertices and not crash b/c of reaching max recursion depth.
    Returns: degSeqGrouped list of pairs (node, target degree)
    '''
    ## degSeq = list of (node, degree)
    degSeqGrouped = []
    group = []
    ##currentTargetDeg
    cTD = degSeq[0][1] ## deg of largest deg node
    for i in range(len(degSeq)):
        if len(group) < k:        ## if group too small
            group.append(degSeq[i])
        elif i > len(degSeq)-k-1:   ## if there arent enough to make another full group, add to prev group
            group.append(degSeq[i])
        elif cTD - degSeq[i][1] > degSeq[i][1] - degSeq[i+k][1]: ## Cmerge > Cnew -> create new group
            degSeqGrouped.append(group)
            group = [degSeq[i]]
            cTD = degSeq[i][1]
        else: 
            group.append(degSeq[i])
    degSeqGrouped.append(group)
    return degSeqGrouped ## [(node,maxdeg)]

def pseudoNodeAnon(G:nx.Graph, k=None, doGreedy=False):
    '''
    create recursion on deg sequence, minimise max deficiency
    incremental alg that goes left to right (assumes 1 indexing) on the deg sequence
    adds the next xth rightmost degree and get best k-partitioning. if there are < 2k degs in sequence, then only 1 partition
    Returns: Ganon anonymous nx.Graph
    '''
    if k == None:
        k = 10
    Ganon = G.copy()
    d = dict(G.degree())
    degSeq = sorted(d.items(), key = lambda item : item[1], reverse=True)  ## list of (node, degree)
    if len(degSeq) < 2*k:
        degDef = {node[0]:degSeq[0][1] - node[1] for node in degSeq}
        md = degSeq[0][1] - degSeq[-1][1] ## difference between largest and smallest degrees
    elif doGreedy:
        degSeqGrouped = greedyGrouping(degSeq, k)
        degDef = {partition[j][0]:partition[0][1]-partition[j][1] for partition in degSeqGrouped for j in range(len(partition))}
        md = max(degDef.values()) ## maximum deficiency
    else:
        degSeqGrouped = []
        i = G.number_of_nodes()
        degSeq = [0] + degSeq
        cost = [None for i in degSeq] ## list of length n, with None as values
        start = [None for i in degSeq]## list of length n, with None as values
        while i>1:
            starti = Start_unspec(degSeq,start,cost,k,i)
            degSeqGrouped.append(degSeq[starti:i+1])    ## need to add 1, b/c we want to include i
            i = starti - 1
        degDef = {partition[j][0]:partition[0][1]-partition[j][1] for partition in degSeqGrouped for j in range(len(partition))}
        ## ^ {node: targetDeg - nodeDeg} 
        degSeq = degSeq[1:]
        md = max(degDef.values()) ## maximum deficiency
        ##  partition is [(node,maxdeg),...,(node,mindeg))] so degDef is nodej: maxdeg-degj
    maxNodeID = max(G.nodes)
    newNodes = list(range(maxNodeID+1, maxNodeID+1 + (1+max(md,k)%2 + max(md,k))))
    Ganon.add_nodes_from(newNodes) ## add [n,...,n+m] to [0,...,n-1]
    ## iterate through new nodes, and add edge to old nodes -> cycle through these slower
    newNodeID = 0
    for (node, deg) in degSeq:
        for i in range(degDef[node]):
            Ganon.add_edge(node, newNodes[newNodeID%len(newNodes)])
            newNodeID +=1
    d = int(newNodeID/len(newNodes)) + 1
    targetDeg = [partition[0][1] for partition in degSeqGrouped]
    if (d in targetDeg and d-1 in targetDeg) or newNodeID == 0: ## also k-anon if all have same deg
        return Ganon
    ## vertices with deg d-1 pair them and add an edge between them
    leftoverNodes = newNodes[newNodeID%len(newNodes):]      ## list of length m-td mod m with nodes of degree d-1
    random.shuffle(leftoverNodes)
    for i in range(0,len(leftoverNodes)-1,2):
        Ganon.add_edge(leftoverNodes[i],leftoverNodes[i+1])
    ## that is [newNodeID to end of list] of new nodes -> if m-td mod m is even and m >= k then we are done
    if len(leftoverNodes)%2 ==0:
        return Ganon
    ## if m-td mod d is odd, then we leave out 1 vertex -> m-1 even and >=k 
    ## add an edge from r to two other new vertices (from full list) -> have degree d+1
    edgesDeg_d = newNodes[0:newNodeID%len(newNodes)] + leftoverNodes[:-1]   ## all but the last element
    random.shuffle(edgesDeg_d)
    Ganon.add_edges_from([(leftoverNodes[-1],edgesDeg_d[0]),(leftoverNodes[-1],edgesDeg_d[1])])
    for i in range(2,len(edgesDeg_d)-1,2):  ## pair off, but els 0 and 1 are already taken
        Ganon.add_edge(edgesDeg_d[i],edgesDeg_d[i+1])
    return Ganon