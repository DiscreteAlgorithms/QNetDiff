from array import array
from audioop import reverse
from cmath import isclose, log10
from code import interact
from distutils.file_util import write_file
from posixpath import split
import string
from sys import stderr
from turtle import color
import pandas as pd

def Calc_Jaccard_SymDiff( bacterias, stageToEdgeToWeight ):
    stages = list( stageToEdgeToWeight.keys() )
    
    bactToJaccard = dict()
    bactToSymDiff = dict()
    for b in bacterias:
        edgeweight_setsum = 0.
        edgeweight_setproduct = 0.
        X = stages[0]
        Y = stages[1]
        for another_b in bacterias:
            # 枝のタプルは常に辞書順で (小さいもの、大きいもの) の順で表す
            b1 = min( b, another_b )
            b2 = max( b, another_b )
            if b1==b2:
                continue
            x_weight = stageToEdgeToWeight[X][b1,b2] if (b1,b2) in stageToEdgeToWeight[X] else 0
            y_weight = stageToEdgeToWeight[Y][b1,b2] if (b1,b2) in stageToEdgeToWeight[Y] else 0

            edgeweight_setsum += max( x_weight, y_weight )
            edgeweight_setproduct += min( x_weight, y_weight )
        if edgeweight_setsum == 0:
            print(f"error: {b}'s setsum = 0")
            exit()
        bactToJaccard[b] = edgeweight_setproduct / edgeweight_setsum
        bactToSymDiff[b] = edgeweight_setsum - edgeweight_setproduct
    return bactToJaccard, bactToSymDiff

def Calc_NetShift_withIsolated( bacterias, stageToEdgeToWeight):
    control, case = list( stageToEdgeToWeight.keys() )
    ctrl_edges = stageToEdgeToWeight[control].keys()
    case_edges = stageToEdgeToWeight[case].keys()

    deg = {b:0 for b in bacterias}
    for v1, v2 in case_edges:
        deg[v1] += 1
        deg[v2] += 1
    max_degree_in_case = max( deg.values() )
    bactToNESH_score = dict()
    for v1 in bacterias:
        intersect = 0.0
        union = 0.0
        unique_in_case = 0.0
        for v2 in bacterias:
            if (v1,v2) in ctrl_edges or (v2,v1) in ctrl_edges:
                weight_in_ctrl = 1
            else:
                weight_in_ctrl = 0
            
            if (v1,v2) in case_edges or (v2,v1) in case_edges:
                weight_in_case = 1
            else:
                weight_in_case = 0

            intersect += min( weight_in_ctrl, weight_in_case ) 
            union += max( weight_in_ctrl, weight_in_case )
            if weight_in_ctrl == 0 :
                unique_in_case += weight_in_case
        X = intersect / union
        Y = unique_in_case / max_degree_in_case
        Z = unique_in_case / union
        bactToNESH_score[v1] = 1 - X + Y + Z
    return bactToNESH_score
    
def Calc_ContinuousNetShift( bacterias, stageToEdgeToWeight ):
    control, case = list( stageToEdgeToWeight.keys() )
    ctrl_edges = stageToEdgeToWeight[control].keys()
    case_edges = stageToEdgeToWeight[case].keys()

    deg = {b:0 for b in bacterias}
    for v1, v2 in case_edges:
        weight = stageToEdgeToWeight[case][ v1, v2 ]
        deg[v1] += weight
        deg[v2] += weight
    max_degree_in_case = max( deg.values() )

    
    bactToNESH_score = dict()
    for v1 in bacterias:
        common = 0.0
        union = 0.0
        unique_in_case = 0.0
        for v2 in bacterias:
            if (v1,v2) in ctrl_edges:
                weight_in_ctrl = stageToEdgeToWeight[control][ v1, v2 ]
            elif (v2,v1) in ctrl_edges:
                weight_in_ctrl = stageToEdgeToWeight[control][ v2, v1 ]
            else:
                weight_in_ctrl = 0
            
            if (v1,v2) in case_edges:
                weight_in_case = stageToEdgeToWeight[case][ v1, v2 ]
            elif (v2,v1) in case_edges:
                weight_in_case = stageToEdgeToWeight[case][ v2, v1 ]
            else:
                weight_in_case = 0

            common += min( weight_in_ctrl, weight_in_case ) 
            union += max( weight_in_ctrl, weight_in_case )
            if weight_in_ctrl == 0 :
                unique_in_case += weight_in_case
        X = 1 - common / union
        Y = unique_in_case / max_degree_in_case
        Z = unique_in_case / union
        bactToNESH_score[v1] = X+Y+Z
    return bactToNESH_score

def Step7( MWU_RANK_TH, CONTRACT_SKIP, PICKUP_SKIP, G, focus_stage, p_val, bacteria_focused, abundance_mean, stageToBactToAbdavg ):
    stages = list( G.keys() )
    bacterias, _= G[stages[0]]

    s0 = stages[0]
    s1 = stages[1]

    degree = dict()
    unique_degree = dict()
    for s in stages:
        degree[s] = dict()
        unique_degree[s] = dict()
        _, E = G[s]
        another_stage = s1 if s == s0 else s0
        _, another_stage_E = G[another_stage]
        for b in bacterias:
            degree[s][b] = 0
            unique_degree[s][b] = 0    
            for another_b in bacterias:
                b1 = min(b,another_b)
                b2 = max(b,another_b)    
                if (b1,b2) in E :
                    degree[s][b] += E[b1, b2] 
                    if (b1,b2) not in another_stage_E :
                        unique_degree[s][b] += E[b1, b2]
    
    non_focus_stage = s0 if focus_stage == s1 else s1

    # degree_diff = { b:degree[focus_stage][b] - degree[non_focus_stage][b] for b in bacterias }
    jaccard, symm_diff = Calc_Jaccard_SymDiff( bacterias, stageToEdgeToWeight={ s0:G[s0][1], s1:G[s1][1] } )

    weighted_nesh = Calc_ContinuousNetShift( bacterias, stageToEdgeToWeight={ s0:G[s0][1], s1:G[s1][1] } )
    nesh = Calc_NetShift_withIsolated( bacterias, stageToEdgeToWeight={ s0:G[s0][1], s1:G[s1][1] } )

    if PICKUP_SKIP:
        if CONTRACT_SKIP:
            print("error: Don't skip both contract and pickup", file=stderr)

    return symm_diff, nesh, weighted_nesh, degree


#keyToThicksColor,
def Output_MultiBar( filepath, keys, param_labels, label_to_param_v, xtick_list=None ):
    BAR_WIDTH = 0.25
    bar_left_v = [i+1 for i in range(len(keys)) ]

    import matplotlib.pyplot as plt
    from utility import ColorCode, color_dict

    plt.figure()
    
    for label, colorCode, index in zip(param_labels, list(ColorCode), range(len(param_labels)) ) :
        plt.bar( 
            [ l + index*BAR_WIDTH for l in bar_left_v], 
            label_to_param_v[label], 
            color=color_dict[colorCode][1], 
            edgecolor=color_dict[colorCode][0] , 
            width=BAR_WIDTH-0.01, align='center'
        
        )
    if xtick_list == None:
        plt.xticks( [ l + BAR_WIDTH/2 for l in bar_left_v], list(keys), rotation=90 )
    else:
        plt.xticks( [ l + BAR_WIDTH/2 for l in bar_left_v], xtick_list, rotation=90 )

    plt.subplots_adjust( bottom = 0.35, top=0.99 )
    plt.savefig( filepath )

def Calc_ContinuousDistribution( bacterias, bactToParam ):
    max_val = max(bactToParam.values())
    min_val = min(bactToParam.values())
    
    NUM_OF_PART = 100
    now_l = min_val




