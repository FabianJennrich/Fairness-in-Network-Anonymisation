## Functions for interacting with outputs from ANONET
## Built to use files output from:
## de Jong, Rachel G., Mark P. J. van der Loo, and Frank W. Takes. "The anonymization problem in social networks." arXiv preprint arXiv:2409.16163 (2024). doi: https://doi.org/10.48550/arXiv.2409.16163


## Imports 
import os
import networkx as nx
import pandas as pd

## function for creating ANONET compatible graphs
'''
the graphs need to have the format

!n=[numNodes]
!m=[numEdges
u:v w
v:u 
w:u

where n is the number of nodes, m the number of edges, and the edges (u,v) are represented twice once u:v and once v:u
'''
def createANONETgraph(workingGraph:nx.Graph, fileName:str):
    '''
    writes to fileName an anonet input graph of workingGraph
    write to anonet compatible file
    adj list format but there is a delimiter (:) btwn source and target
    '''
    sorter = dict(zip(sorted(workingGraph.nodes()),
                      range(workingGraph.number_of_nodes())))
    workingGraph = nx.relabel_nodes(workingGraph, sorter)
    with open(fileName, 'w') as f:
        numNodes = "!n=" + str(workingGraph.number_of_nodes()) + "\n"
        f.write(numNodes)
        numNodes = "!m=" + str(workingGraph.number_of_edges()) + "\n"
        f.write(numNodes)
        for node in sorted(workingGraph.nodes()):
            nautyline = str(node) + ":" + " ".join([str(ns) for ns in sorted([n for n in workingGraph.neighbors(node)])]) + ";\n"
            f.write(nautyline)


## and for reading ANONET output graph and partitions
def read_anonetoutput(GTgraph:nx.Graph, path:str):
    '''read path into networkx graph'''
    ## recreates original node names for anonet graph
    desorter = dict(zip(range(GTgraph.number_of_nodes()),
                        sorted(GTgraph.nodes())))
    ## remove all : and ; tokens to get regular adj list format
    with open(path, 'r') as f:
        adjList = f.readlines()
    adjList = [line.translate(str.maketrans("","",":;")) for line in adjList]
    anonGraph = nx.parse_adjlist(adjList, nodetype=int)
    return nx.relabel_nodes(anonGraph, desorter)

def read_anonet_partition(GTgraph:nx.Graph, path:str):
    '''
    Assumes that each line in the input file is a partition, each entry in the line is a node ID. 
    Returns: partitionDict {node: partitionID} dictionary of which partition each node is in,
        kanonDict {node: kanonValue} dictionary of the size of the partition each node is in
    '''
    ## recreates original node names for anonet graph
    desorter = dict(zip(range(GTgraph.number_of_nodes()),
                        sorted(GTgraph.nodes())))
    with open(path,'r') as f:
        partitionList = f.readlines()
    partitionList = [line.translate(str.maketrans("","",",\n")).split() for line in partitionList]
    partitionDict = {node: partitionID for partitionID in range(len(partitionList)) for node in partitionList[partitionID]}
    kanonDict = {node:len(partition) for partition in partitionList for node in partition}
    partitionDict = {desorter[node]: partitionDict[node] for node in partitionDict}
    kanonDict = {desorter[node]: kanonDict[node] for node in kanonDict}
    return partitionDict, kanonDict