import networkx as nx

from community import community_louvain



def GetNxGraph( bacterias, adjacency_matrix ):
    stages = list( adjacency_matrix.keys() )
    graph = nx.Graph()
    for b1 in bacterias:
        added = False
        for b2 in bacterias:
            if b1 >= b2 :
                continue
            if (b1,b2) in adjacency_matrix:
                graph.add_edge( b1, b2, weight=adjacency_matrix[b1,b2] )
                added = True
        if not added:
            graph.add_node(b1)
    return graph

def Step5_1( bacterias, stageToAdjacency_matrix):
    stages = list( stageToAdjacency_matrix.keys() )
    stageToNxGraph = dict()
    for s in stages:
        stageToNxGraph[s] = GetNxGraph( bacterias, stageToAdjacency_matrix[s])
    stageToCluster = dict()
    for s in stages:
        stageToCluster[s] = community_louvain.best_partition( stageToNxGraph[s], random_state = 0 )
    return stageToCluster


import union_find 
def Step5_2( bacterias, sup_category, abundance_mean, stageToBactToCluster):
    s1,s2 = tuple( stageToBactToCluster.keys() )
    uf  = union_find.UnionFind(bacterias)
    for b1 in bacterias:
        for b2 in bacterias:
            if b1 >= b2:
                continue
            if uf.isSame( b1, b2 ):
                continue
            term1 = stageToBactToCluster[s1][b1] == stageToBactToCluster[s1][b2] and stageToBactToCluster[s2][b1] == stageToBactToCluster[s2][b2]
            term2 = sup_category[b1] == sup_category[b2]
            if term1 and term2:
                uf.unite( b1, b2 )

    parentToMost = dict()
    for b in bacterias:
        parent = uf.find(b)
        if parent not in parentToMost.keys():
            parentToMost[ parent ] = b
        else:
            parentToMost[ parent ] = max( parentToMost[parent], b, key = lambda x: abundance_mean[x] )

    contract_set_list = dict()
    for b in bacterias:
        most = parentToMost[uf.find(b)]
        if most not in contract_set_list.keys():
            contract_set_list[most] = []
        contract_set_list[most].append(b)
    return contract_set_list


def Step5_3( contract_set_list, stageToAdjacency_matrix ):
    bacteria_contracted = contract_set_list.keys()
    adjacency_matrix_contracted = dict()
    for s in stageToAdjacency_matrix.keys():
        adjacency_matrix_contracted[s] = dict()
        for b1 in bacteria_contracted:
            for b2 in bacteria_contracted:
                if b1 >= b2:
                    continue
                contract_set1 = contract_set_list[ b1 ]
                contract_set2 = contract_set_list[ b2 ]
                weight_sum = 0. 
                cnt = 0
                for v1 in contract_set1:
                    for v2 in contract_set2:
                        if (v1,v2) not in stageToAdjacency_matrix[s]:
                            continue
                        weight_sum += stageToAdjacency_matrix[s][v1,v2]
                        cnt += 1
                if cnt == 0:
                    continue
                adjacency_matrix_contracted[s][ b1, b2 ] = weight_sum / cnt
    return bacteria_contracted, adjacency_matrix_contracted
        