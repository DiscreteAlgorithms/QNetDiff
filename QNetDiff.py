from array import array
from cmath import isclose
from os import path
from scipy import stats
import sys
import math

import utility

from read_input import ReadCntTable_Whole, ReadCategory
from step0 import Step0_CalcRelativeAbundance
from step3 import Step3
from For_SparCC import Run_SparCC
from step1 import Step1
from step2 import Step2_1, Step2_2, Step2_3
from step4 import Get_Induced_InBothGroups, GetNeighbors, Step4
from step5 import Step5

# mann-whitney-u(MWU) p_val threshold 
# bacteria whiche's p_val is lower than this TH is regarded "significantly" increased
MWU_P_TH = 0.005

# MWU rank threshold 
# how many bacterias will be picked up by rank of MWU p-val
MWU_RANK_TH = 10

EDGE_THRESHOLD = 0.4

output_folda_path = "QNetDiff_output/"

# tableFile :file of relative abundance table. size is [number of bacteria] x [numbar of sample]
# sampleFile:file of sample infomation each has [count of sample] line

def CalcAllAbundanceMean( abundance ):
    bacterias = list( abundance.keys() )
    abundance_mean = dict()
    for b in bacterias:
        sum = 0.
        cnt = 0 
        for s in abundance[b].keys():
            for abd in abundance[b][s]:
                cnt += 1
                sum += abd
        if cnt == 0:
            abundance_mean[ b ] = 0
        else:
            abundance_mean[ b ] = sum / cnt
    return abundance_mean

def Output_ContractInfo( parentToContractSet, bact_to_sup_category ): 
    infos = []
    parents = list(parentToContractSet.keys())
    contract_group_num = 0
    for parent in parents:
        s = ""
        contract_set = parentToContractSet[parent]
        if len(contract_set) == 1:
            continue
        contract_group_num += 1
        for b in contract_set:
            s += f"{b}, "
        s += f"\n are contracted to {parent} (sup_category:{bact_to_sup_category[parent]}, group size: {len(contract_set)})\n"
        infos.append(s)
    infos = sorted(infos)
    with open("QNetDiff_output/contractinfo.txt","w")as f:
        f.write(f"number of representative bacteria: {len(parents)}\n")
        f.write(f"number of similar group (size >= 2): {contract_group_num}\n")
        for s in infos:
            f.write(s)

def AbundanceToCircleSize( abundance, min_abundance, circlesize_min, grad  ):
    from math import log10
    return  circlesize_min + grad*( log10(abundance/min_abundance).real )

from nx_utlility import Show_net
from utility import similar_color_dict, ColorCode
def Output_NxGraph( filename_footer, groupToVE_Tuple, to_sort_feature, bacteria_focused, bact_to_abundance_mean, draw_label=True, node_size_change=True, edge_width_change=True, sort_vertices=True, edgeweight_range=[0.0,1,0] ):
    groups = groupToVE_Tuple.keys()
    edge_count = dict()
    for s in groups:
        _, E = groupToVE_Tuple[s]
        for e in E:
            if e in edge_count:
                edge_count[e] += 1
            else:
                edge_count[e] = 1 
    
    import networkx as nx
    NODEGROUP_FOCUSED = 1
    NODEGROUP_UNFOCUSED = 0
    
    for s in groups:
        V, E = groupToVE_Tuple[s]
        if sort_vertices:
            sorted_v = sorted( V, key=lambda v:to_sort_feature[v], reverse=True )
        else:
            sorted_v = V
        g = nx.Graph()
        for v in sorted_v:
            if v in bacteria_focused:
                g.add_node( v, group=NODEGROUP_FOCUSED )
            else:
                g.add_node( v, group=NODEGROUP_UNFOCUSED )
        for e in E:
            g.add_edge( e[0], e[1], weight = E[e] )
        
        swapped_edge_count = dict()
        for e in g.edges:
            if e in edge_count:
                swapped_edge_count[e] = edge_count[e]
            else:
                swapped_edge_count[e] = edge_count[(e[1],e[0])]
        min_abundance = min([bact_to_abundance_mean[b] for b in V])

        if node_size_change:
            size_vec = [
                    AbundanceToCircleSize( bact_to_abundance_mean[v], min_abundance, 40, 350 )
                    for v in V
                ]
        else:
            size_vec = [ 20 for v in V]
        Show_net( 
            g, 
            output_folda_path+f"{s}_net{filename_footer}.png",
            node_to_size = size_vec,
            nodeGroup_to_color = { 
                NODEGROUP_FOCUSED:   ColorCode.RED.value,
                NODEGROUP_UNFOCUSED: ColorCode.GRAY.value
            }, 
            # nodeGroup_to_framecolor = {
            #     NODEGROUP_FOCUSED:   ColorCode.RED.value,
            #     NODEGROUP_UNFOCUSED: ColorCode.GRAY.value
            # },
            edgecolor_v = [ 
                ColorCode.GREEN.value if swapped_edge_count[(u,v)] == 1 else ColorCode.GRAY.value
                for u, v in g.edges 
            ],
            draw_label=draw_label,
            edge_width_change=edge_width_change,
            edge_weight_min=edgeweight_range[0], edge_weight_max=edgeweight_range[1]
        )

def Output_GraphInfo( path, bacterias, groupToEdgeToWeight ):
    groups = list( groupToEdgeToWeight )
    with open( path, "w" ) as f:
        f.write(f"vertex_num:{len(bacterias)}\n")
        for s in groups:
            weight_sum = 0.0
            edge_num = 0
            f.write(f"======= in {s} ======\n")
            for b1 in bacterias:
                for b2 in bacterias:
                    if b1 >= b2 :
                        continue
                    weight = 0
                    if (b1,b2) in groupToEdgeToWeight[s]:
                        weight = groupToEdgeToWeight[s][b1,b2]
                    elif (b2,b1) in groupToEdgeToWeight[s]:
                        weight = groupToEdgeToWeight[s][b2,b1]
                    
                    weight_sum += weight
                    if weight > EDGE_THRESHOLD:
                        edge_num += 1
            f.write(f"edge_num:{edge_num}\n")
            f.write(f"edge_weight_sum:{weight_sum}\n")

def GetParentAndContractV( bacteria, contract_set_list ):
    for parent in contract_set_list:
        contract_set = contract_set_list[parent]
        if bacteria in contract_set:
            return parent, contract_set
    return None

FIG_LINE_WIDTH = 0.3

def Convert_Feature_to_Size( feature_dict, size_range ):
    max_val = max( list(feature_dict.values()) )
    min_val = min( list(feature_dict.values()) )

    grad = (size_range[1]-size_range[0]) / (max_val-min_val)
    return_dict = dict()
    for key, val in feature_dict.items():
        return_dict[key] = size_range[0] + (val-min_val)*grad
    return return_dict

def GetNx_g( V, edgeToWeight ):
    import networkx as nx
    nx_g = nx.Graph()
    nx_g.add_nodes_from( V )
    for e in edgeToWeight:
        nx_g.add_edge( e[0], e[1], weight = edgeToWeight[e] )
    return nx_g

def main():
    
    args = sys.argv
    
    if len(args) < 5 or 7 < len(args):
        print("arg size must be 5, 6 or 7 ", file=sys.stderr)
        return  

    for i in range(2):
        if not path.isfile(args[1]):
            v = ["table","sample"]
            print(f"{utility.ToOrdinal(i)} arg is not file path (it must be {v[i]} file path)", file=sys.stderr )
            return
    
    tableFilePath = args[1]
    sampleFilePath = args[2]
    categoryFilePath = args[3]
    focusGroup = args[4]

    corr_file_name_without_group_and_extension = ""
    if len(args) >= 6 :
       corr_file_name_without_group_and_extension = args[5] + "_"

    SPARCC_SKIP = False
    if len(args) == 7 :
        if args[6] == 'SKIP_SPARCC':
            SPARCC_SKIP = True
        else:
            print( "If you want to skip making correlation file, final arg must be \"SKIP_SPARCC\". ", file=sys.stderr )
            exit()

    import pandas as pd
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("display.precision", 6)

    bacterias, groups,  bactToGroupToCnts, groupToSamples = ReadCntTable_Whole( tableFilePath, sampleFilePath )
    
    bactToGroupToAbds = Step0_CalcRelativeAbundance( bacterias, groups, bactToGroupToCnts )
    print("step0 finished")
    
    if args[4] not in groups:
        print(f"Third arg:{args[4]} is not found from group column in sampleFile.", file=sys.stderr)
        return

    if groups[0]==focusGroup:
        baseGroup = groups[1]
    else:
        baseGroup = groups[0]

    if SPARCC_SKIP:
        print("sparCC skipped")
    else:
        Run_SparCC( groups, bacterias, bactToGroupToCnts, groupToSamples )
    
    groupToEdgeToWeight, bacteria_effective, groupToVTupleToCorr = Step1( EDGE_THRESHOLD, groups, bacterias)
    print("step1 finished")

    abundance_mean = CalcAllAbundanceMean( bactToGroupToAbds )

    sup_category = ReadCategory( categoryFilePath )

    groupToBactToCluster = Step2_1( bacteria_effective, groupToEdgeToWeight )

    parentToContractSet = Step2_2( bacteria_effective, sup_category, abundance_mean, groupToBactToCluster )

    bacteria_contracted, groupToEdgeToWeight_contracted = Step2_3( parentToContractSet, groupToEdgeToWeight )
    print("step2 finished")

    Output_ContractInfo( parentToContractSet, sup_category )

    bactToGroupToAbds = { b:bactToGroupToAbds[b] for b in bacteria_contracted }
    bacteria_focused, p_val = Step3( MWU_P_TH, MWU_RANK_TH, baseGroup, focusGroup, bactToGroupToAbds )
    print("step3 finished")

    groupTo_VE_tuple = Step4( bacteria_focused, bacteria_contracted, groupToEdgeToWeight_contracted )
    print("step4 finished")

    V = groupTo_VE_tuple[groups[0]][0]
    groupToEdgeToWeight_picked = {s:groupTo_VE_tuple[s][1] for s in groupTo_VE_tuple.keys()}

    abundance_mean_inGraph = {b:abundance_mean[b] for b in V}

    groupToBactToAbdAvg = { s:{ b: sum(bactToGroupToAbds[b][s]) / len(bactToGroupToAbds[b][s]) for b in V } for s in groups}
    sym_diff, groupToBactToDegree = Step5( groupTo_VE_tuple )
    print("step5 finished")

    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = 'Bahnschrift'
    plt.rcParams["font.size"] = 10
    
    # グラフ出力
    Output_NxGraph( "", groupTo_VE_tuple, sym_diff, bacteria_focused, abundance_mean_inGraph )

    # 各指標の表出力
    df_dict = {
        "QNetDiff": sym_diff,
        "p-value":{b:p_val[b] for b in V},
        f"degree_in_{groups[0]}":groupToBactToDegree[groups[0]],
        f"degree_in_{groups[1]}":groupToBactToDegree[groups[1]],
        f"abundance_average_in_{groups[0]}":groupToBactToAbdAvg[groups[0]],
        f"abundance_average_in_{groups[1]}":groupToBactToAbdAvg[groups[1]]
    }
    output_df = pd.DataFrame(df_dict)
    print( output_df.sort_values("QNetDiff", ascending=False) , file=open(output_folda_path+f"{groups[1]}-{groups[0]}_result_table.txt","w") )


if __name__ == "__main__":
    main()