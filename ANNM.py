## A New Noise Node Method
## Implemented based on the pseudo code in:
## S. Hamzehzadeh and S. M. Mazinani, ‘‘ANNM: A new method for adding noise nodes which are used recently in anonymization methods in social networks,’’ Wireless Pers. Commun., vol. 107, no. 4, pp. 1995–2017, Aug. 2019.

## ANNManon(G:nx.Graph, k = 10)
## produces k=10 degree anonymity
## Line number comments refer to lines in the pseudo code of the ANNM paper

import numpy as np
import networkx as nx
from itertools import groupby

def degreeGroupsCreation(G:nx.Graph):
    '''
    groups the degree sequence of the graph G by grouping nodes with the same degree value
    Input: G nx.Graph
    Returns: degSeq = [(node,deg)] sorted highest deg first,
        degSeqGrouped = [[(node,deg1),...,(node,deg1)],...,[(node,degn),...,(node,degn)]] sorted by descending degree
    '''
    d = dict(G.degree())
    degSeq = sorted(d.items(), key = lambda item : item[1], reverse=True)  ## list of (node, degree)
    degSeqGroupsed = [[pair for pair in value] for key,value in groupby(degSeq, lambda x: x[1])]    ## lists of (node,deg) each list has only 1 deg
    return degSeq, degSeqGroupsed

def degreePrioritization(G:nx.Graph, degSeqGrouped, k):
    '''
    Input: G nx.Graph to be anonymized,
        degSeqGrouped the node degree pairs as list sorted highest degree first of lists grouped by degree value [[(node,deg1),...,(node,deg1)],...,[(node,degn),...,(node,degn)]]
    Returns: insecureGroupPriorities = {groupID: priority} where groupID is the index of the group in degSeqGrouped
    '''
    insecureGroupPriorities = {}
    for i in range(len(degSeqGrouped)):
        s = len(degSeqGrouped[i])
        if s <= k-1:
            insecureGroupPriorities[i] = int((k-s <= s/2))   
            ## evaluates to 1 if true, 0 if D > S/2
    return insecureGroupPriorities

def targetGroupsCreation(G:nx.Graph, P, degSeqGrouped, k):
    '''
    Groups insecure nodes such that the degree anonymity is larger than k. Using members of groups with priority 0, we complete groups priority 1. If there are less than k insecure nodes, add all insecure nodes to smallest secure group.
    Input: G nx.Graph,
        P is the priority of insecure groups
        degSeqGrouped grouped degree sequence
        k the minimum group size in the anonymous graph
    Returns: Groups_prime dictionary of groupID: [(node, deg),...(node,deg)], 
        gPO a list of nodes ordered by priority (1 then 0)
    '''
    ## P = {groupID: priority}
    ## iterate through them backwards, so that they are sorted lowest deg to highest
    groupsPriority1 = [key for key,value in P.items() if value == 1][::-1]
    groupsPriority0 = [key for key,value in P.items() if value == 0][::-1]
    gPO = groupsPriority1 + groupsPriority0 ## list reordering nodes

    Groups   = {i:degSeqGrouped[gPO[i]] for i in range(len(gPO))}   ## newID : [(node, deg)]
    G_p = {i:P[gPO[i]] for i in range(len(gPO))}               ## newID : priority
    Groups_prime = {}

    allInsecNodes = [nodePair for group in Groups.values() for nodePair in group]
    if len(allInsecNodes) < k:
        minSecGroup = sorted([group for group in degSeqGrouped if len(group)>=k], key= len)[0]
        return {0:allInsecNodes + minSecGroup}, 0

    order = len(gPO)
    existGroups0 = len(groupsPriority0)
    for i in range(order):
        if G_p[i] == 1 and existGroups0:
            Groups_prime[i] = Groups[i]     ## 19
            Groups[i] = [] ## empty it so that we can keep easy track of non-added elements later
            ## starting at i+1 guarantees j>i
            for j in range(i+1,order):      
                if G_p[j] == 0:
                    ## moves on to next j if the first one is empty
                    while len(Groups_prime[i]) <=k and Groups[j] != []:   
                        Groups_prime[i].append(Groups[j][0])
                        if len(Groups[j]) == 1:
                            existGroups0 -= 1
                        Groups[j] = Groups[j][1:]
        elif G_p[i] == 1 and not existGroups0:
            Groups_prime[i] = Groups[i]
            Groups[i] = []
            for j in range(i+1, order): ## for j to order, if (j>i)
                while len(Groups_prime[i]) <= k and Groups[j] != []:
                    Groups_prime[i].append(Groups[j][0])
                    Groups[j] = Groups[j][1:]
        else: ## if P[gPO[i]] == 0:
            Groups_prime[i] = Groups[i]
            for j in range(i+1, len(gPO)):
                while len(Groups_prime[i]) <= k and Groups[j] != []:
                    Groups_prime[i].append(Groups[j][0])
                    Groups[j] = Groups[j][1:]       
    leftoverNodes = [pair for group in Groups.values() for pair in group]
    if len(leftoverNodes)>k-1:
        Groups_prime[order] = leftoverNodes
    else:
        ## add to maximal group with priority 1
        lastIndexp1 = max([i for i in range(len(gPO)) if P[gPO[i]]==1])
        Groups_prime[lastIndexp1] += leftoverNodes
    return Groups_prime, gPO

def addingNoiseAnon(G:nx.Graph, group_prime):
    '''
    Adds noise nodes <<to compensate for the difference in degree between the group degree and nodes added in the third step>>
    Input: G nx.Graph in process of anonymization
        group_prime dictionary of target groupID: nodes list
    Returns: Ganon nx.Graph anonymized graph
    '''
    bc = dict(nx.betweenness_centrality(G))
    bcset = sorted(bc.items(), key = lambda item : item[1])  ## list of (node, bc) lowest to highest
    bco = {bcset[i][0]:i for i in range(len(bcset))}      ## vertex pos in bcset, so node: pos
    gp_working = {groupID:nodes for (groupID,nodes) in group_prime.items() if nodes != []}
    group_prime_bcSorted = {k: sorted(v, key = lambda item: bc[item[0]]) for (k,v) in gp_working.items()}
    ## should be a dict of same form as group_prime, except sorted by bc, groupID : [(node,deg)]
    
    noisetag = True
    noisecounter = 0        ## 4
    Ganon = G.copy()
    newNodes = {} ## = {degree:[noise nodes w/ that target degree]}

    ## Adding noise nodes and anonymization
    for groupID in group_prime_bcSorted:    ## for each group:
        ## consider node with lowest betweenness centrality -> iterate through my sorted list
        ## ^ 22, 26, 34 
        group = group_prime_bcSorted[groupID]
        targetDeg = max([deg for (node,deg) in group])
        if targetDeg not in newNodes:
            newNodes[targetDeg] = []
        for (node,deg) in group: ## 23 for(i=0 to group'.s)
            ## if group'.V[i].deg < group'.deg & noisetag == if deg < targetDeg & noisetag
            if noisetag:      ## 24
                for j in range(targetDeg - deg): ## 27 for (j=0 to targetDeg-deg)
                    ## ^ completes the check deg < targetDeg
                    newNodeID = Ganon.number_of_nodes()
                    newNodes[targetDeg].append(newNodeID)
                    Ganon.add_node(newNodeID)
                    Ganon.add_edge(newNodeID, node) ## connect to the node
                    noisecounter += 1
                    noisetag = False
            else:            ## 32
                for new_node in sorted(newNodes[targetDeg], key = lambda x: Ganon.degree(x)):
                    ## add to nodes with lowest deg first
                    if Ganon.degree(node) >= targetDeg:     ## do not add edge if unneccessary
                        continue
                    if Ganon.degree(new_node) < targetDeg:
                        ## maximum degree of noise nodes is target degree   -> make sure there is not a massive degree difference
                        Ganon.add_edge(new_node, node)
                for j in range(targetDeg - Ganon.degree(node)):     ## if existing nodes are insufficient, add new ones
                    ## for(j=0 to group'deg-Vi.deg)
                    ## if that is not enough noise nodes add more from node with next lowest betweenness centrality
                    newNodeID = Ganon.number_of_nodes()
                    newNodes[targetDeg].append(newNodeID)
                    Ganon.add_node(newNodeID)
                    Ganon.add_edge(newNodeID, node)
                    noisetag = True

    ## try to connect such that all nodes are k-anon
    ## if not possible, move to node with next lowest betweenness centrality and continue
    ## kanon noise nodes
    possibleDegs = set(dict(Ganon.degree()).values())
    unsecNoiseNodes = {} ## {nodeID: degDef}
    for deg in newNodes:
        for node in newNodes[deg]:
            if Ganon.degree(node) not in possibleDegs:
                ## ^ If the degree of created noise node [not] exists 
                ## among the degrees are available in degree groups of social network graph
                unsecNoiseNodes[node] = deg - Ganon.degree(node)    ## node: deg deficiency

    if len(unsecNoiseNodes) == 0:
        return Ganon

    unsecNodesSorted = sorted(unsecNoiseNodes.items(),key = lambda item: item [1], reverse=True)        ## sorted by deg def    
    targetDegSeqUnsec = list(unsecNoiseNodes.values())
    mycounter = 0
    while not nx.is_graphical(targetDegSeqUnsec, method="hh") and mycounter <=100:
        try:
            min_odd_Target = sorted([deg for deg in newNodes if deg%2 == 1])[0]
            targetDegSeqUnsec.append(min_odd_Target)
            mycounter+=1
        except:
            mycounter = 200
    if mycounter >= 100:
        for i in range(len(unsecNodesSorted)):
            nodei = unsecNodesSorted[i][0]
            for (nodej,ddefj) in unsecNodesSorted[i:]:
                if unsecNoiseNodes[nodei] >0 and unsecNoiseNodes[nodej] >0:
                    Ganon.add_edge(nodej, nodei)
                    unsecNoiseNodes[nodej] -= 1
                    unsecNoiseNodes[nodei] -= 1
            unsecNodesSorted = sorted(unsecNoiseNodes.items(),key = lambda item: item [1], reverse=True)    ## resort -> still checks 0 entries
        ## does not actually guarantee k-anon for noise nodes
        return Ganon
    ## if it is successful in creating realizable deg seq:
    connGraph = nx.havel_hakimi_graph(targetDegSeqUnsec)
    connGraphNodesSorted = sorted(dict(connGraph.degree()).items(), key= lambda item:item[1], reverse=True)
    trans = {}
    moreNoiseNodes = []
    for i in range(len(connGraphNodesSorted)):
        (oldID, oldDeg) = unsecNodesSorted[i-len(moreNoiseNodes)]
        (newID, newDeg) = connGraphNodesSorted[i]
        if newDeg == oldDeg:
            trans[newID] = oldID
        elif newDeg > oldDeg:
            plusID = Ganon.number_of_nodes()+1+len(moreNoiseNodes)
            moreNoiseNodes.append(plusID)
            trans[newID] = plusID
    Ganon.add_nodes_from(moreNoiseNodes)
    for (u,v) in connGraph.edges():
        Ganon.add_edge(trans[u], trans[v])
    return Ganon

def ANNManon(G:nx.Graph, k = 10):
    '''
    Creates the anonymous graph, with 10-degree anonymity. 
    Input: G nx.Graph to be anonymized,
        k the minimum degree anonymity
    Returns: Ganon anonymous nx.Graph
    '''
    degSeq, degSeqGrouped = degreeGroupsCreation(G)
    insecGroupP = degreePrioritization(G, degSeqGrouped, k)
    group_prime, gPO = targetGroupsCreation(G,insecGroupP,degSeqGrouped, k)
    Ganon = addingNoiseAnon(G, group_prime)
    return Ganon