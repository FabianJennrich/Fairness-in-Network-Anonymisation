## Pipeline for Graph Anonymization and Evaluation
## Load graphs to anonymize
## Anonymize
## Save anonymous graphs to file
## 

## Imports
import pandas as pd ## for the sensible organisation of privacy score outcomes
import networkx as nx   ## used in all the anon methods and metrics
import datetime     ## used to time anonymizations

## Get the anonymization algorithms

from ANNM        import ANNManon   ## k=10
from DPFT        import DP_FT
from edgeEnt     import sample_uncertain_graph
from ENDP        import projectionOrderedEdgeInsertion     ## threshold 80% quantile deg
from FiM         import fim_anon ## 50 fim nodes
from KDVEM       import KDVEMapproach     ## k = 10
from linkMirage  import linkMirage   ## t = 5 -> length of random walk
from linkPerturb import transform  ## t = 5 -> same as linkMirage
from OLP         import optimizationLinkPerturbation  ## beta = 0.1
from PBCN        import PBCNanonymisation  ## t = 10 -> num groups, K = edges/100 -> num noise edges
from pseudoNode  import pseudoNodeAnon   ## k = 10
from randomWalk  import rwAnon       ## t = 5 -> length of random walk -> same param vals as linkMirage
from subgraphDP  import perturbGraph ## k = 5 -> subgraph depth, beta = 0.5 -> ratio of edges that are fake
from twoStep     import probing ## k = 10
from upperApprox import PPGPworkingEx
from kDegSeq     import k_degSeq_Anon ## k = 10, max failure count (in a row) = 20

## Import evaluation metrics
from myMetrics import graph_evaluation

## code for running and saving the graph anonymization
def testAnon(G:nx.Graph):
    '''
    This is a test function that returns the input graph unchanged.
    Useful for checking that the anonMethodDict is working
    Also useful to check that none of the anonymisation methods changed the ground truth graph.
    '''
    return G

## dictionary of anonymization methods
anonMethodDict = {'ANNM': ANNManon,
                  'DPFT':DP_FT,
                  'edgeEnt':sample_uncertain_graph,
                  'ENDP': projectionOrderedEdgeInsertion,
                  'FiM': fim_anon,
                  'KDVEM': KDVEMapproach,
                  'linkMirage': linkMirage,
                  'linkPerturb': transform,
                  'OLP':optimizationLinkPerturbation,
                  'PBCN': PBCNanonymisation,
                  'pseudoNode': pseudoNodeAnon,
                  'randomWalk': rwAnon,
                  'subgraphDP': perturbGraph,
                  'upperApprox': PPGPworkingEx,
                  'twoStep': probing,
                  'kDegSeq': k_degSeq_Anon,
                  'Naive' : testAnon}

def graph_anonymisation(G:nx.Graph, graphName:str, methodDict = anonMethodDict, date:str = "NoDateGiven"):
    '''
    Anonymizes the input graph G according to the methods given in methodDict
    saves them as adjacency list in anonymousGraphs/graphName_method_date
    Returns: anonGraphsDict dictionary of form {method: Ganon}
    '''
    ## dictionary of anon algorithms with a unique name string for the key
    ## and the function as the value
    anonGraphsDict = {}
    for method in methodDict.keys():
        anonGraph_path = 'anonymousGraphs/' + graphName + '_' + method + "_" + date 
        ##try:    ## if a method does not work simply move on to the next one
        start = datetime.datetime.now()
        anonGraph = methodDict[method](G)
        ## save anonGraph to file
        nx.write_adjlist(anonGraph, anonGraph_path)
        anonGraphsDict[method] = anonGraph
        runtime = datetime.datetime.now() - start
        print(runtime, anonGraph_path)
        ##except:
        ##    print(anonGraph_path, 'Did not generate.')
    return anonGraphsDict

## function for running and saving graph evaluation metrics
def saveGraphEval(GTgraph:nx.Graph, AnonGraphDict:dict, graphName:str, metricNames = [], date:str = "NoDateGiven"):
    '''
    run all metrics in myMetrics on each of the anon graphs, 
    save the results as csv file metricsResults/graphName_[metricNames]_date.csv
    '''
    myDF = graph_evaluation(GTgraph, AnonGraphDict, metricNames)
    filepath = 'metricsResults/' + graphName + "_" + str(metricNames) + "_" + date + '.csv'
    myDF.to_csv(filepath)


## loading in the anonymous graphs from file, assumes that the ANONET graphs were generated seperately and placed in the same folder.
anonMethods = ['ANNM', 'ANONET-deg','ANONET-count',
               'ANONET-isom','DPFT','edgeEnt',
               'ENDP','FiM', 'KDVEM', 
               'linkMirage', 'linkPerturb', 'OLP', 
               'PBCN', 'pseudoNode', 'randomWalk', 
               'subgraphDP', 'upperApprox', 'twoStep',
               'kDegSeq', 'Naive']

def getAnonymousGraphs(graphName:str, date:str, methods = anonMethods):
    '''
    Loads in anonymous graphs from adjacency list with nodetype int. The file location for each anonymity method is assumed to be "anonymousGraphs/Gname_method_date".
    Returns anonGraphsDict {method: anonGraph} dictionary of anonymous graphs
    '''
    anonGraphsDict = {}
    for method in methods:
        anonGraph_path = 'anonymousGraphs/' + Gname + '_' + method + "_" + date
        try:
            anonGraph = nx.read_adjlist(anonGraph_path, nodetype = int)
            anonGraphsDict[method] = anonGraph
            print(Gname, method," was successfully imported")
        except:
            print(anonGraph_path, 'was missing')
    return anonGraphsDict

## Code for running in commandline
## while loop that runs for user input
## input graph, graph name, date (to identify different anon of same graph)

while True:
    if str(input("Anonymize a Graph (y/n):")).strip() in ['y','Y','yes','Yes']:
        Gpath = str(input("Path to graph to anonymize (must be gml with integer nodes):"))
        G = nx.read_gml(Gpath, destringizer = int)
        Gname = str(input("Name of the Graph:"))
        date = str(input("Date:"))

        anonGraphs = graph_anonymisation(G, 
                                         graphName= Gname, 
                                         methodDict= anonMethodDict, 
                                         date= date)
    elif str(input("Evaluate Anonymous Graphs (y/n):")).strip() in ['y','Y','yes','Yes']:
        Gpath = str(input("Path to ground truth graph (must be gml with integer nodes):"))
        G = nx.read_gml(Gpath, destringizer = int)
        Gname = str(input("Name of the Graph:"))
        date = str(input("Date:"))
        anonGraphsDict = getAnonymousGraphs(Gname, date)
        saveGraphEval(G, anonGraphsDict, Gname, date = date)
    elif str(input("Quit (y/n):")).strip() in ['y','Y','yes','Yes']:
        break