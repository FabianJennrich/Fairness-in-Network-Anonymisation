## myMetrics
## defining graph metrics for analyzing the anonymous graphs

## imports
import networkx as nx
import pandas as pd
import numpy as np
from collections import Counter


def k_degree_anonymity(G:nx.Graph, GanonDict:dict):
    '''
    Computes the k-degree anonymity of each node for each of the graphs in GanonDict.
    Input: G the non-anonymous graph (unused, present for consistency),
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs,
    Returns: resultsDict {graphName: {node: k-degree anon}} dictionary of dictionaries of k-degree anonymity for each node for each anonymous graph
    '''
    ## sort pairs by value, groupby value, get lengths of groups
    resultsDict = {}
    for graphID in GanonDict:
        k_dict = dict(GanonDict[graphID].degree())   ## to adapt this to k-(something)-anon, change this dictionary
        anonDict = {}                       ## dictionary that contains the size of each equivalence class
        for node in GanonDict[graphID].nodes():              ## for each node
            equivClass = k_dict[node]       ## get the equiv class (degree) of the node
            if equivClass in anonDict:      ## check if the equivalence class exists in my dictionary keys
                anonDict[equivClass] = anonDict[equivClass] + 1     ## if yes, increment
            else: 
                anonDict[equivClass] = 1    ## if not add it
        resultsDictID = graphID + " k-Deg Anon"
        resultsDict[resultsDictID] = {keyNode:anonDict[k_dict[keyNode]] for keyNode in GanonDict[graphID].nodes()}
    return resultsDict

def jointDegSequence(G:nx.Graph): 
    '''returns list of edge endpoint degrees'''
    return [(min(G.degree(u), G.degree(v)), max(G.degree(u), G.degree(v))) for (u,v) in G.edges() ]

def jointDegSequence_subset(G:nx.Graph, edges, removing=None): 
    '''returns list of edge endpoint degrees'''
    d = dict(G.degree())
    if removing != None:
        for u in G.neighbors(removing):
            d[u] -= 1
    return [(min(d[u],d[v]), max(d[u],d[v])) for (u,v) in edges]

def distribution_difference_from_counter(dist1, dist2, returnDict = False): 
    '''
    returns the difference between two distributions given as dictionaries
    value between 0(same dict) and sum(dict1)+sum(dict2) (no keys in common)
    '''
    ## start computation of distance
    mydistance = {key: dist1[key]-dist2[key] for key in set(dist1).union(dist2)}
    if returnDict:
        return mydistance
    return  sum(np.abs(list(mydistance.values())))

def differential_privacy(G:nx.Graph, GanonDict:dict):
    '''
    Computes the differential privacy (DP) score of each of the graphs in GanonDict.
    Input: G the non-anonymous graph, 
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs, 
    Returns: resultsDict {graphName: {node: DP score}} dictionary of dictionaries of DP score for each node for each anonymous graph
    '''
    GT_jdDist = Counter(jointDegSequence(G))
    d = G.degree()
    resultsDict = {} ## output dictionary
    GT_jdDist_etorem = {}
    GT_jdDist_etoadd = {}
    for node in G.nodes():
        edges_N1 = G.edges(G.neighbors(node))
        edges_N1minNode = [e for e in edges_N1 if e not in G.edges(node)]
        GT_jdDist_etorem[node] = Counter(jointDegSequence_subset(G, edges_N1)) ## edges to remove
        GT_jdDist_etoadd[node] = Counter(jointDegSequence_subset(G,edges_N1minNode,removing=node)) ## edges to add
    for graphID in GanonDict:
        resultsDictName = graphID + " Diff Priv"
        ## get the distribution for the anon graph
        Anon_jdDist = Counter(jointDegSequence(GanonDict[graphID]))
        GTminAnon_dict = distribution_difference_from_counter(GT_jdDist,Anon_jdDist, returnDict=True)
        ## where GT is the positive values, and Anon is the (-) values
        dG1_S = sum(np.abs(list(GTminAnon_dict.values())))
        ## dont bother with other computation if the anonymised graph is the 'same' as the GT graph
        if dG1_S == 0:  
            resultsDict[resultsDictName] = {key:np.inf for key in G.nodes()}
            continue
        dG2_S = {}
        for node in G.nodes():
            helper = GTminAnon_dict.copy() ## since we have GT-anon, we can add and subtract from the dict ...
            for key in GT_jdDist_etorem[node]:
                helper[key] -= GT_jdDist_etorem[node][key]
            for key in GT_jdDist_etoadd[node]:
                if key in helper:
                    helper[key] += GT_jdDist_etoadd[node][key]
                else:
                    helper[key] = GT_jdDist_etoadd[node][key]
            dG2_S[node] = sum(np.abs(list(helper.values()))) ## ... and then sum at the end
            
        diff_priv_dict = {node: np.inf if dG2_S[node] == 0 else max(dG1_S/dG2_S[node], dG2_S[node]/dG1_S) for node in G.nodes()}
        resultsDict[resultsDictName] = diff_priv_dict
    return resultsDict

def lcc(G:nx.Graph, GanonDict:dict):
    '''
    Computes the local clustering coefficient (LCC) for each of the nodes of each of the graphs in GanonDict using nx.clustering.
    Input: G the non-anonymous graph (unused, present for consistency), 
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs, 
    Returns: resultsDict {graphName: {node: LCC(node)}} dictionary of dictionaries of LCC for each node for each anonymous graph
    '''
    resultsDict = {}
    for graphID in GanonDict:
        resultsDictNames = graphID + " Local Clustering Coeff."
        lcc = nx.clustering(GanonDict[graphID])
        resultsDict[resultsDictNames] = lcc
    return resultsDict

def get_bc(G:nx.Graph, GanonDict:dict):
    '''
    Computes the betweenness centrality (BC) for each of the nodes of each of the graphs in GanonDict using nx.betweenness_centrality, with log(n+1) pivots, where n is the number of nodes in the anonymous graph.
    Input: G the non-anonymous graph (unused, present for consistency), 
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs, 
    Returns: resultsDict {graphName: {node: BC(node)}} dictionary of dictionaries of BC for each node for each anonymous graph
    '''
    resultsDict = {} 
    for graphID in GanonDict:
        numPivots = int(np.log(GanonDict[graphID].number_of_nodes()+1)+ 1) 
        resultsDictName = graphID + " Betweenness Centrality"
        if numPivots == 1:
            bc = dict(nx.betweenness_centrality(GanonDict[graphID]))
        else: bc = dict(nx.betweenness_centrality(GanonDict[graphID],k=numPivots))
        resultsDict[resultsDictName] = bc
    return resultsDict

def get_degree(G:nx.Graph, GanonDict:dict):
    '''
    Computes the degree for each of the nodes of each of the graphs in GanonDict using nx.Graph.degree
    Input: G the non-anonymous graph (unused, present for consistency), 
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs, 
    Returns: resultsDict {graphName: {node: deg(node)}} dictionary of dictionaries of BC for each node for each anonymous graph
    '''
    resultsDict = {}
    for graphID in GanonDict:            
        resultsDictNames = graphID + " Degree"
        deg = dict(GanonDict[graphID].degree())
        resultsDict[resultsDictNames] = deg
    return resultsDict

def changed_edges(G:nx.Graph, GanonDict:dict):
    '''
    Computes the the number of edges changed during anonymization by summing for each node the number of neighbors it has gained and the number of neighbors it has lost
    Input: G the non-anonymous graph, 
        GanonDict {graphName: nx.Graph} dictionary of anonymous graphs, 
    Returns: resultsDict {graphName: {node: deg(node)}} dictionary of dictionaries of BC for each node for each anonymous graph
    '''
    ## works by getting {node: [neighbors]} for each graph, and comparing the elements
    N1GT = {n:set(G.neighbors(n)) for n in G.nodes()}

    resultsDict = {}
    for graphID in GanonDict:
        resultsDictName = graphID + " Changed Edges"
        N1AG = {n: set(GanonDict[graphID].neighbors(n)) if n in GanonDict[graphID].nodes() else set([]) for n in G.nodes()}
        diffAG = {n: len(N1GT[n].union(N1AG[n])) - len(N1GT[n].intersection(N1AG[n])) for n in N1GT}
        ## use GT to iterate, so that we consider the relevant nodes
        resultsDict[resultsDictName] = diffAG
    return resultsDict

def graph_evaluation(G:nx.Graph, AnonGraphDict:dict, metricNames:list = [], graphName:str = "MetricsResults"):
    '''
    Inputs: G Ground truth graph with community labels, 
        AnonGraphDict {graphName: nx.Graph} dictionary of anonymised graphs, 
        metricNames list of the metrics to be used as strings,
        graphName string for identifying G in filename
    Returns: total_df pandas dataframe with results of several graph anonymity metrics.
        Saves total_df to file as "metricsResults/graphName[metricNames].csv"
    '''
    ## requires that all metrics take the same inputs in the same order
    ## that is metric(GT_graph, AnonGraphsDict)
    metricsDict = {'kdeg': k_degree_anonymity, 
                   'DP': differential_privacy,
                   'lcc': lcc,
                   'bc': get_bc,
                   'deg':get_degree,
                   'edge_d':changed_edges,
                  }
    
    if metricNames == []:   ## if we enter an empty list run all of the metrics
        metricNames = list(metricsDict.keys())

    commValuebyNode = nx.get_node_attributes(G, 'community')
    commDict = {}
    for node in commValuebyNode:
        if commValuebyNode[node] in commDict:
            commDict[commValuebyNode[node]].append(node)
        else:
            commDict[commValuebyNode[node]] = [node]

    outputDFs = []
    for metric in metricNames:
        outputDFs.append(pd.DataFrame.from_dict(metricsDict[metric](G, AnonGraphDict)))
        print(metric, " concluded")

    total_df = pd.concat(outputDFs, axis = 1)
    try: ## attempt to add community value labels
        total_df = pd.concat([total_df,
                              pd.DataFrame.from_dict({"Community Value": dict(G.nodes(data='community'))})],
                             axis = 1)
    except:
        print("Failed to add community labels")
        
    total_df.to_csv("metricsResults/"+graphName+str(metricNames)+".csv")
    return total_df