## Friend in the Middle
## Implemented based on the description in:
## F. Beato, M. Conti, and B. Preneel, “Friend in the middle (FiM): Tackling de-anonymization in social networks,” in Proc. 5th Int. Workshop Security Soc. Netw., San Diego, CA, USA, 2013, pp. 279–284.

## fim_anon(G:nx.Graph)

import networkx as nx
import random

def fim_anon(G:nx.Graph):
    ''' Anonymizes G using Friend in the Middle process. Chooses 50 nodes from G to be FiM nodes, and then for each node reroutes 0.2 of its edges through one of the FiM nodes. 
    Returns: Ganon anonymous nx.Graph'''
    fim_nodes = random.sample(list(G.nodes()),k=50)
    ## can construct it from largest to smallest deg node, i.e.
    ## can iterate through the FiM nodes in order, so that we dont get duplicates
    ## we need percentage of max deg FiM nodes
    ## for each node's neighbors, we add either that edge, or an FiM edge -> O(n^2) or O(n*h) h avg deg
    edgeList = []
    nonEdgeList = []
    fim_iterator = 0
    for n in G.nodes():
        for v in G.neighbors(n):
            if random.random() <= 0.2:
                edgeList.append((n,fim_nodes[fim_iterator]))
                edgeList.append((fim_nodes[fim_iterator],v))
                nonEdgeList.append((n,v))
            else:
                edgeList.append((n,v))
            fim_iterator = (fim_iterator + 1)%50
    Ganon = nx.Graph()
    Ganon.add_nodes_from(G.nodes())
    Ganon.add_edges_from(edgeList)
    Ganon.remove_edges_from(nonEdgeList)
    return Ganon