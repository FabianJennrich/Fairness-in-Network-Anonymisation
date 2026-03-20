## Random Walk based Anonymisation
## Implemented based on the pseudo code in:
## Y. Guo, Z. Liu, Y. Zeng, R. Wang, and J. Ma, ‘‘Preserving privacy for hubs and links in social networks,’’ in Proc. Int. Conf. Netw. Netw. Appl. (NaNA), Oct. 2018, pp. 263–269.

## rwAnon(G:nx.Graph)
## max random walk length of 5

import numpy as np
import networkx as nx
import random

def randomWalk(G:nx.Graph, u, distance):
    '''Returns the terminal vertex of a random walk of length distance starting at u'''
    vertex = u
    for i in range(distance):
        vertex = random.choice(list(G.neighbors(vertex)))
    return vertex

def degreePerturbation(G:nx.Graph, tc, Mc=10):
    G1 = G.copy()
    nthresh = np.quantile(sorted([d for n,d in G1.degree()], reverse=True),0.1)    ## threshold for top 10% of degrees
    N = [n for n,d in G1.degree() if d >= nthresh]   ## nodes with degree in top 10%
    for u in N:
        count = 1
        while True:
            v = randomWalk(G1,u,tc)
            degu = G1.degree(u)
            degv = G1.degree(v)
            count += 1
            if count > Mc:       ## continue if count <= M
                break
            ## continue if degv > degu
            if degv <= degu and not [node for node in G1.neighbors(u) if node != v and node not in G1.neighbors(v)] == []:    
                break
        if count <= Mc:
            Nv = G1.neighbors(v)
            Nucond = [node for node in G1.neighbors(u) if node != v and node not in Nv]
            ## Nu = randomly choose (degu - degv) nodes from neighbors of u s.t.
            ## v not in Nu and Nu intersection Nv = 0
            sampleSize = min(degu-degv, len(Nucond))
            Nu = random.sample(Nucond, k=sampleSize)
            for z in Nu:
                G1.remove_edge(u,z)
                G1.add_edge(v,z)
    return G1

def linkPerturbation(G1:nx.Graph,t,M0=10):
    Gprime = nx.Graph()
    for u in G1.nodes():
        number = 1
        for v in G1.neighbors(u):
            count = 1
            while True:
                z = randomWalk(G1,v, t-1)
                count += 1
                if count > M0:  ## continue if count small enough
                    break
                if u != z and (u,z) not in Gprime.edges():  ## continue if u=z or (u,z) already exists
                    break
            if count <= M0:
                if number == 1:
                    Gprime.add_edge(u,z)
                else:
                    prob = (0.5 * G1.degree(u) - 1) / (G1.degree(u) - 1)
                    if random.random() <= prob:
                        Gprime.add_edge(u,z)
        number += 1
    return Gprime

def rwAnon(G:nx.Graph):
    ## tc and t are random walk lengths
    tc = 5
    t = 5
    G1 = degreePerturbation(G, tc) 
    return linkPerturbation(G1, t)