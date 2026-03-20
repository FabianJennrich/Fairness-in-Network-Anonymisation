## Link Perturbation Anonymisation
## Implemented based on the pseudo code in:
## P. Mittal, C. Papamanthou, and D. Song, “Preserving link privacy in social network based systems,” in Proc. 20th Annu. Netw. Distrib. Syst. Security Symp. (NDSS), San Diego, CA, USA, 2013, pp. 1–15.

## transform(G:nx.Graph, t=5, M=10)
## random walk length t=5, and number of attempts per edge M=10

import networkx as nx
import random

def randomWalk(G:nx.Graph, u, distance):
    '''Returns the terminal vertex of a random walk of length distance starting at u'''
    vertex = u
    for i in range(distance):
        vertex = random.choice(list(G.neighbors(vertex)))
    return vertex

def transform(G:nx.Graph, t=5, M=10):
    '''Perturb undirected graph G using perturbation t and maximum loop count M
    t is the random walk length -> chose 5 from paper results
    M is the number of tries to find a suitable edge'''
    Gprime = nx.Graph()
    Gprime.add_nodes_from(G.nodes())    ## vertices stay the same
    for u in G.nodes(): ## foreach u in G
        count = 1
        for v in G.neighbors(u):    ## foreach neighbor v of u
            loop = 1
            ## do ... until (u=z OR (u,z) in G') AND (loop <=M)
            while True:
                z = randomWalk(G, v, t-1)   ## t-1 hop rand walk from v
                loop += 1
                if loop > M:    ## stop if do loop too often
                    break
                if u != z and (u,z) not in Gprime.edges():  ## this should be correct
                    break
            if loop <= M:
                if G.degree(u) > 1:
                    prob = (0.5 * G.degree(u) - 1) / (G.degree(u) - 1)
                else:
                    prob = 0
                if count == 1:
                    Gprime.add_edge(u,z)
                elif random.random() <= prob:
                    Gprime.add_edge(u,z)
                count += 1      ## increase count only if an edge was added -> else move on to the next neighbor without increasing count
    return Gprime