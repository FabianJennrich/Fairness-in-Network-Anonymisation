## Edge iNsertion based Differential Privacy
## Implemented based on the pseudo code in:
## K. R. Macwan and S. J. Patel, ‘‘Node differential privacy in social graph degree publishing,’’ Procedia Comput. Sci., vol. 143, pp. 786–793, Jan. 2018.

## projectionOrderedEdgeInsertion(G:nx.Graph)

import numpy as np
import networkx as nx

def projectionOrderedEdgeInsertion(G:nx.Graph):
    ''' Anonymizes G via ENDP. The graph is reconstructed by imposing a degree maximum. The 80th percentile of degrees is taken as the degree max.
    Returns: Ganon anonymous nx.Graph'''
    theta = np.quantile([d for n,d in G.degree()], 0.8)
    Ganon = nx.Graph()
    Ganon.add_nodes_from(G.nodes())
    E_theta = []
    d = {node: 0 for node in G.nodes()}
    Vprime = sorted(G.nodes(), key = lambda node: G.degree(node), reverse=False)   
    ## sorted node list based on deg [(node,maxdeg),...,(node,mindeg)]
    for v in Vprime:
        if d[v] <theta:
            neigh_v = G.neighbors(v) 
            for u in sorted(neigh_v, key = lambda node:d[node]):
                if d[u] < theta and d[v] < theta and not Ganon.has_edge(u,v):
                    Ganon.add_edge(u,v)
                    d[u] += 1
                    d[v] += 1
    return Ganon