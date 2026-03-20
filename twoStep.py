## Two Step Anonymisation
## Implemented based on the pseudo code in:
## K. Liu and E. Terzi, “Towards identity anonymization on graphs,” in Proc. ACM SIGMOD Int. Conf. Manag. Data, Vancouver, BC, Canada, 2008, pp. 93–106.

## probing(G:nx.Graph, k = 10, useGreedy = False)
## Note: partitioning is done recursively and gives a max recursion error in large graphs. In graphs where this is an issue, useGreedy uses greedy partitioning of the degree sequence

import numpy as np
import networkx as nx
import random

def funcI(degSeq, i, j):
    ''' returns the result of function I(d[i,j]) defined in the paper. We assume degSeq of form [(nodeID, nodeDeg)] sorted highest deg to lowest deg. If i and j are indices of deqSeq then:
    I(d[i,j]) = sum_{l=i}^{j} (d(i)-d(l)) -> do j+1 so that index j is included
    '''
    return degSeq[i][1]*(j+1-i) - sum([item[1] for item in degSeq[i:j+1]])

def greedyGrouping(degSeq, k):
    '''a greedy method for determining the partition of degSeq
    this should allow it to run with more vertices and not crash b/c of reaching max recursion depth.
    '''
    ## degSeq = list of (node, degree)
    degSeqGrouped = []
    degDict = {}
    group = []
    ##currentTargetDeg
    cTD = degSeq[0][1] ## deg of largest deg node
    for i in range(len(degSeq)):
        if len(group) < k:        ## if group too small
            group.append(degSeq[i])
            degDict[degSeq[i][0]] = cTD
        elif i +1+k > len(degSeq)-1:   ## if there arent enough to make another full group, add to prev group
            group.append(degSeq[i])
            degDict[degSeq[i][0]] = cTD
        else:
            Cmerge = (cTD - degSeq[i][1]) + funcI(degSeq, i+1, i+1+k)
            Cnew = funcI(degSeq, i, i+k)
            if Cmerge > Cnew: ## make new group starting at i
                degSeqGrouped.append(group)
                group = [degSeq[i]]
                cTD = degSeq[i][1]
                degDict[degSeq[i][0]] = cTD
            else: 
                group.append(degSeq[i])
                degDict[degSeq[i][0]] = cTD
    degSeqGrouped.append(group)
    return degDict

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
    return degSeq[x][1]*(y-x) - sum([i[1] for i in degSeq[x:y+1]])

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

def degreeAnon(d, k:int):
    '''take a degree sequence, and return a degree sequence such that the anon cost is minimized
    d = [(nodeID, deg)] such that d[1]>= ... >= d[n]
    Returns a dictionary of nodeID: targetDegree'''
    ## Da(d[1:i]) = deg anon cost of d[1:i] and I(d[i:j]) = deg anon cost when all i, i+1, ..., j are in the same anon group
    ## I(d[i,j]) = sum l = i to j (d(i)-d(l))
    ## for i <= 2k Da(d[1,i]) = I
    ## for i >= 2k Da(d[1:i]) = min { min over k<=t<= i-k {Da}}
    if len(d) < 2*k:
        dhat = {d[i][0]:d[-1][1] for i in d.keys()}
        return dhat
    i = len(d)
    d = [0] + d ## add element for 1 indexing
    DA = [None for node in d] ## costs
    start = [None for node in d]
    degSeqGrouped = []
    while i>1:
        starti = Start_unspec(d,start,DA,k,i)
        degSeqGrouped.append(d[starti:i+1])    ## need to add 1, b/c we want to include i
        i = starti - 1
    d_hat = {partition[j][0]: partition[0][1] for partition in degSeqGrouped for j in range(len(partition))}
    d = d[1:] ## remove element for 1 indexing
    return d_hat

def priority(d, Ggt:nx.Graph):
    '''
    Allows the construction of deg anon graphs with similar high edge intersection directly 
    without using greedy_swap.
    Similar to constructGraph, but makes two passes over sorted deg seq.
    First pass, considers only nodes such that (v,vprime) in E.
    If not enough, does second pass such that (v,vprime) not in E'''
    G = nx.Graph()
    G.add_nodes_from(d.keys())  ## get the node set
    if sum([deg for deg in d.values()]) %2 ==1:  ## if sum of degrees is odd
        print("failed, sum odd")
        return False, None
    while True:
        if any([deg<0 for deg in d.values()]):
            print([(k,v) for k,v in d.items() if v<0])
            return False, None
        if all([deg==0 for deg in d.values()]):
            return True, G
        v = random.choice([n for n,deg in d.items() if deg > 0])    ## random node with d(v) >0
        ## vprime s.t. (v,vprime) in E
        vprime = sorted([(n,d[n]) for n in Ggt.neighbors(v) if (d[n]>0 and n!=v)], key = lambda x:x[1], reverse = True) 
        if len(vprime) < d[v]:
            vprime += sorted([(n,deg) for n,deg in d.items() if (deg>0 and (n not in Ggt.neighbors(v)) and n!=v)], key= lambda x:x[1], reverse=True)
        Vdv = vprime[:d[v]]
        d[v] = 0
        for (w,dw) in Vdv:
            if w == v: print(v)
            G.add_edge(v,w)
            d[w] -= 1

def probing(G:nx.Graph, k = 10, useGreedy = False):
    '''Anonymises G'''
    d = sorted(dict(G.degree()).items(), key = lambda item : item [1], reverse= True)
    if useGreedy: dhat = greedyGrouping(d,k)
    else: dhat = degreeAnon(d,k) ## returns dict node:targetDeg
    n = G.number_of_nodes()
    realizable, Ghat = priority(dhat, G)
    while not realizable:
        noise = random.choice(range(1, 5))
        d[-1] = (d[-1][0], min(d[-1][1] + noise, n-1)) ## must replace the tuple in its entirety
        ## repartition and get target deg sequence
        d = sorted(d, key = lambda item : item [1], reverse= True)
        if useGreedy: dhat = greedyGrouping(d,k)
        else: dhat = degreeAnon(d,k)
        realizable, Ghat = priority(dhat, G)
    return Ghat