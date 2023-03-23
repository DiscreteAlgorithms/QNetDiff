from sqlite3 import connect


def GetNeighbors( bacteria_focused, whole_bacterias, stageToEdges ):
    ret = []
    
    for b in whole_bacterias:
        if b in bacteria_focused:
            ret.append( b )
            continue
        connected = False
        for s in stageToEdges.keys():
            edges = stageToEdges[s]
            for another_b in bacteria_focused:
                if (b, another_b) in edges or (another_b, b) in edges:
                    connected = True
                    break
            if connected:
                break
        if connected:
            ret.append(b)
    return ret

def Get_Induced_InBothStages( bacterias, stageToEdgeToWeight ):
    ret_stageToEdgeToWeight = dict()
    for s in stageToEdgeToWeight.keys():
        ret_stageToEdgeToWeight[s] = dict()
        for b1 in bacterias:
            for b2 in bacterias:
                if (b1, b2) not in stageToEdgeToWeight[s]:
                    continue
                ret_stageToEdgeToWeight[s][ b1, b2 ] = stageToEdgeToWeight[s][ b1,b2 ]
    return ret_stageToEdgeToWeight


def Step6( PICKUP_SKIP, bacteria_focused, whole_bacterias, stageToEdgeToWeight ):
    V = []
    stages = list(stageToEdgeToWeight.keys())
    stageToEdges = { s:stageToEdgeToWeight[s].keys() for s in stages }
    if PICKUP_SKIP:
        V = list(whole_bacterias)
    else:
        V = GetNeighbors( bacteria_focused, whole_bacterias, stageToEdges )

    G = dict()
    for s in stages:
        E = dict()
        for v1 in V:
            for v2 in V:
                if (v1,v2) not in stageToEdges[s]:
                    continue
                weight = stageToEdgeToWeight[s][ v1, v2 ]
                if weight > 0 :
                    E[ v1, v2 ] = weight
        G[s] = (V,E)
    return G
