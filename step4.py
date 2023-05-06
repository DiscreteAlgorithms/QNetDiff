def GetNeighbors( bacteria_focused, whole_bacterias, groupToEdges ):
    ret = []
    
    for b in whole_bacterias:
        if b in bacteria_focused:
            ret.append( b )
            continue
        connected = False
        for s in groupToEdges.keys():
            edges = groupToEdges[s]
            for another_b in bacteria_focused:
                if (b, another_b) in edges or (another_b, b) in edges:
                    connected = True
                    break
            if connected:
                break
        if connected:
            ret.append(b)
    return ret

def Get_Induced_InBothGroups( bacterias, groupToEdgeToWeight ):
    ret_groupToEdgeToWeight = dict()
    for s in groupToEdgeToWeight.keys():
        ret_groupToEdgeToWeight[s] = dict()
        for b1 in bacterias:
            for b2 in bacterias:
                if (b1, b2) not in groupToEdgeToWeight[s]:
                    continue
                ret_groupToEdgeToWeight[s][ b1, b2 ] = groupToEdgeToWeight[s][ b1,b2 ]
    return ret_groupToEdgeToWeight


def Step4( bacteria_focused, whole_bacterias, groupToEdgeToWeight ):
    V = []
    groups = list(groupToEdgeToWeight.keys())
    groupToEdges = { s:groupToEdgeToWeight[s].keys() for s in groups }
    V = GetNeighbors( bacteria_focused, whole_bacterias, groupToEdges )

    G = dict()
    for s in groups:
        E = dict()
        for v1 in V:
            for v2 in V:
                if (v1,v2) not in groupToEdges[s]:
                    continue
                weight = groupToEdgeToWeight[s][ v1, v2 ]
                if weight > 0 :
                    E[ v1, v2 ] = weight
        G[s] = (V,E)
    return G
