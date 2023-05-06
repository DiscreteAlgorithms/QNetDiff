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

def Calc_SymDiff( bacterias, stageToEdgeToWeight ):
    stages = list( stageToEdgeToWeight.keys() )
    
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
        bactToSymDiff[b] = edgeweight_setsum - edgeweight_setproduct
    return bactToSymDiff

def Step5( G ):
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

    symm_diff = Calc_SymDiff( bacterias, stageToEdgeToWeight={ s0:G[s0][1], s1:G[s1][1] } )

    return symm_diff, degree





