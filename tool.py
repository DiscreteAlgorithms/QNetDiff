from array import array
from cmath import isclose
from os import path
from sqlite3 import connect
from turtle import color
from unicodedata import category
from scipy import stats
import sys
import math

import utility

from read_input import ReadCntTable_Whole, ReadCategory
from step1 import Step1
from step2 import Step2
from step3 import Step3
from step4 import Step4
from step5 import Step5_1, Step5_2, Step5_3
from step6 import Get_Induced_InBothStages, GetNeighbors, Step6
from step7 import Step7

# mann-whitney-u(MWU) p_val threshold 
# bacteria whiche's p_val is lower than this TH is regarded "significantly" increased
MWU_P_TH = 0.005

# MWU rank threshold 
# how many bacterias will be picked up by rank of MWU p-val
MWU_RANK_TH = 10

EDGE_THRESHOLD = 0.4

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

def Output_ContractInfo( contract_set_list, bact_to_sup_category ): 
    infos = []
    contracted_num = 0
    contract_group_num = 0
    for parent in contract_set_list.keys():
        s = ""
        contract_set = contract_set_list[parent]
        if len(contract_set) == 1:
            continue
        contracted_num += len(contract_set)-1
        contract_group_num += 1
        for b in contract_set:
            s += f"{b}, "
        s += f"\n are contracted to {parent} (sup_category:{bact_to_sup_category[parent]}, group size: {len(contract_set)})\n"
        infos.append(s)
    infos = sorted(infos)
    with open("tool_output/contractinfo.txt","w")as f:
        f.write(f"contracted bacteria(vertex) num: {contracted_num}\n")
        f.write(f"contract group num: {contract_group_num}\n")
        for s in infos:
            f.write(s)

def Output_ForCytoscape( foldapath, fileNameHeader, V, stageToEdgeToWeight, bacteria_focused, node_attirbute_dict=None ):
    s0, s1 = stageToEdgeToWeight.keys()
    
    for s in [s0,s1]:
        another_stage = s0 if s==s1 else s1
        with open( foldapath + f"{fileNameHeader}_{s}.tsv","w") as f:
            f.write("v1\tv2\tweight\tisUnique\n")
            for (v1,v2) in stageToEdgeToWeight[s]:
                f.write(f"{v1}\t{v2}\t{stageToEdgeToWeight[s][v1,v2]}\t{(v1,v2) not in stageToEdgeToWeight[another_stage]}\n")
            for v in V[1:]:
                if (v,V[0]) not in stageToEdgeToWeight[s]:
                    f.write(f"{v}\t{V[0]}\t0.0\tTrue\n")

    with open( foldapath + f"{fileNameHeader}_nodeAttribute.tsv", "w" ) as f:
        f.write("v\tisfocused")
        if node_attirbute_dict != None:
            for attribute_name in node_attirbute_dict.keys():
                f.write(f"\t{attribute_name}")
        f.write("\n")
        for b in V:
            f.write(f"{b}\t{ 'true' if b in bacteria_focused else 'false' }")
            if node_attirbute_dict != None:
                for attribute_name in node_attirbute_dict.keys():
                    if b in node_attirbute_dict[attribute_name]:
                        f.write( f"\t{node_attirbute_dict[attribute_name][b]}" )
            f.write("\n")

def Output_ForNetShiftApp( stageToEdgeList ):
    for s in stageToEdgeList.keys():
         with open(f"for_paper_output/netshift/for_netshift_{s}.tsv", "w")as f:
            for e in stageToEdgeList[s]:
                f.write(f"{e[0]}\t{e[1]}\n")

def AbundanceToCircleSize( abundance, min_abundance, circlesize_min, grad  ):
    from math import log10
    return  circlesize_min + grad*( log10(abundance/min_abundance).real )

def For_tex_output(file_path, df_dict, bacterias, bacteria_focused):
    with open(file_path,"w") as f:
        for feature in df_dict.keys():
            f.write(f" & {feature}")
        f.write(" \\\\ \\hline \hline \n")
        for b in reversed( sorted( bacterias, key = lambda b: list(df_dict.values())[0][b] ) ):
            if bacteria_focused is not None and b in bacteria_focused:
                f.write(f"{{\\it \\color{{red}} {b} }}")
            else:
                f.write(f"{{\\it {b} }}")
            for bactToFeature in df_dict.values():
                str = f"{bactToFeature[b]:.3}"
                if 'e' in str: #指数表記をtexの上付き文字に変換
                    header, footer = str.split('e')
                    str = header + ' \\times 10^{' + footer[0]+footer[2] +'}'
                f.write( " & $" + str + "$")

            f.write(f" \\\\ \\hline\n")


from nx_utlility import Show_net
from utility import similar_color_dict, ColorCode
def Output_NxGraph( filename_footer, stageToVE_Tuple, to_sort_feature, bacteria_focused, bact_to_abundance_mean, draw_label=True, node_size_change=True, edge_width_change=True, sort_vertices=True, edgeweight_range=[0.0,1,0] ):
    stages = stageToVE_Tuple.keys()
    edge_count = dict()
    for s in stages:
        _, E = stageToVE_Tuple[s]
        for e in E:
            if e in edge_count:
                edge_count[e] += 1
            else:
                edge_count[e] = 1 
    
    import networkx as nx
    NODEGROUP_FOCUSED = 1
    NODEGROUP_UNFOCUSED = 0
    
    for s in stages:
        V, E = stageToVE_Tuple[s]
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
            f"tool_output/graphs/{s}_net{filename_footer}.svg",
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

def Output_DifferGraph( V,bacteria_focused, stageToEdgeToWeight, nodeToSize ):
    stages = list( stageToEdgeToWeight.keys() )
    edgeToColor = dict()
    import networkx as nx
    nx_g = nx.Graph()

    NODEGROUP_FOCUSED = 1
    NODEGROUP_UNFOCUSED = 0
    for v in V:
        if v in bacteria_focused:
            nx_g.add_node( v, group=NODEGROUP_FOCUSED )
        else:
            nx_g.add_node( v, group=NODEGROUP_UNFOCUSED )
    min_weight = 10
    max_weight = 0.0
    for e, w0 in stageToEdgeToWeight[stages[0]].items():
        w = None
        if e in stageToEdgeToWeight[stages[1]]:
            w1 = stageToEdgeToWeight[stages[1]][e]
            if w0 > w1:
                w = w0-w1
            pass
        else:
            w = w0
        if w is not None:
            nx_g.add_edge( e[0], e[1], weight=w )
            edgeToColor[e] = ColorCode.BLUE.value
            min_weight = min( min_weight, w )
            max_weight = max( max_weight, w )
        
    for e, w1 in stageToEdgeToWeight[stages[1]].items():
        w = None
        if e in stageToEdgeToWeight[stages[0]]:
            w0 = stageToEdgeToWeight[stages[0]][e]
            if w1 > w0:
                w = w1-w0
            pass
        else:
            w = w1
        if w is not None:
            nx_g.add_edge( e[0], e[1], weight=w )
            edgeToColor[e] = ColorCode.BLUE.value
            min_weight = min( min_weight, w )
            max_weight = max( max_weight, w )

    Show_net( 
        nx_g, 
        f"tool_output/graphs/differ_H-S0_net.svg",
        node_to_size = nodeToSize,
        nodeGroup_to_color = { 
            NODEGROUP_FOCUSED:   ColorCode.RED.value,
            NODEGROUP_UNFOCUSED: ColorCode.GRAY.value
        }, 
        edgecolor_v = [ 
            edgeToColor[u,v] if (u,v) in edgeToColor else edgeToColor[v,u]
            for u,v in nx_g.edges 
        ],
        edge_weight_min=0.0, edge_weight_max=1.0,
        node_shape = 'h'
    )

def Output_GraphInfo( path, bacterias, stageToEdgeToWeight ):
    stages = list( stageToEdgeToWeight )
    with open( path, "w" ) as f:
        f.write(f"vertex_num:{len(bacterias)}\n")
        for s in stages:
            weight_sum = 0.0
            edge_num = 0
            f.write(f"======= in {s} ======\n")
            for b1 in bacterias:
                for b2 in bacterias:
                    if b1 >= b2 :
                        continue
                    weight = 0
                    if (b1,b2) in stageToEdgeToWeight[s]:
                        weight = stageToEdgeToWeight[s][b1,b2]
                    elif (b2,b1) in stageToEdgeToWeight[s]:
                        weight = stageToEdgeToWeight[s][b2,b1]
                    
                    weight_sum += weight
                    if weight > EDGE_THRESHOLD:
                        edge_num += 1
            f.write(f"edge_num:{edge_num}\n")
            f.write(f"edge_weight_sum:{weight_sum}\n")

def Output_GraphInfo_forTeX( path, bacterias, stageToEdgeToWeight ):
    stages = list( stageToEdgeToWeight )
    with open( path, "w" ) as f:
        f.write(f"whole vertices:{len(bacterias)}\n")
        f.write(f"ステージ & エッジを持つノードの数  & エッジの数 & エッジ重みの総和 \\\\ \\hline \n")
        for s in stages:
            weight_sum = 0.0
            edge_num = 0
            f.write(f"{s} & ")
            having_edge_vertice_num = 0
            for b1 in bacterias:
                having_edge = False
                for b2 in bacterias:
                    if b1 >= b2 :
                        continue
                    weight = 0
                    if (b1,b2) in stageToEdgeToWeight[s]:
                        weight = stageToEdgeToWeight[s][b1,b2]
                    elif (b2,b1) in stageToEdgeToWeight[s]:
                        weight = stageToEdgeToWeight[s][b2,b1]
                    
                    weight_sum += weight
                    if weight > EDGE_THRESHOLD:
                        edge_num += 1
                        having_edge = True
                if having_edge:
                    having_edge_vertice_num += 1
            f.write(f"{having_edge_vertice_num} & ")
            f.write(f"{edge_num} & {weight_sum} \\\\ \\hline \n")

def GetParentAndContractV( bacteria, contract_set_list ):
    for parent in contract_set_list:
        contract_set = contract_set_list[parent]
        if bacteria in contract_set:
            return parent, contract_set
    return None

def ResearchOneBact_beforeAndAfterContract( research_bact_name, whole_bacterias, parentToContractSet, stageToEdgeToWeight_before, stageToEdgeToWeight_after, stageToBactToCluster ):
    
    #　縮約前のグラフ
    neighbor_bacts_before = GetNeighbors( [research_bact_name], whole_bacterias, stageToEdgeToWeight_before  ) #including research_bact
    in_graph_bacts = set()
    parents = set()
    bactToParent = dict()
    for b in neighbor_bacts_before:
        parent, contract_v = GetParentAndContractV( b, parentToContractSet )
        parents.add( parent )
        for b2 in contract_v:
            if b2 not in in_graph_bacts:
                in_graph_bacts.add( b2 )
                bactToParent[b2] = parent
    in_graph_bacts = list(in_graph_bacts)
    parents = list(parents)

    before_contract_G = Get_Induced_InBothStages( in_graph_bacts, stageToEdgeToWeight_before )
    Output_ForCytoscape( 
        "for_paper_output/", f"{research_bact_name}_beforeContract", 
        in_graph_bacts , before_contract_G, 
        bacteria_focused = [research_bact_name],
        node_attirbute_dict = { "parent":bactToParent, "cluster_H":stageToBactToCluster["Healthy"], "cluster_S0":stageToBactToCluster["Stage_0"] }
    )

    # 縮約後のグラフ
    after_contract_G = Get_Induced_InBothStages( parents, stageToEdgeToWeight_after )    
    Output_ForCytoscape( 
        "for_paper_output/", f"{research_bact_name}_afterContract", 
        parents, after_contract_G, 
        bacteria_focused = [research_bact_name],
        node_attirbute_dict = { "parent":bactToParent, "cluster_H":stageToBactToCluster["Healthy"], "cluster_S0":stageToBactToCluster["Stage_0"] }
    )

FIG_LINE_WIDTH = 0.3

def Output_pValInSomeStages( bacterias, p_val_dicts, MWU_P_TH ):
    stage_list = list( p_val_dicts.keys() )
    import matplotlib.pyplot as plt
    from math import log10
    _, ax = plt.subplots()
    plots = []
    bacterias = sorted( bacterias, reverse=True, key=lambda b:-log10(p_val_dicts[stage_list[-1]][b]).real )
    bacterias_plot = []
    for b in bacterias:
        if min( [p_val_dicts[s][b] for s in stage_list] ) < MWU_P_TH:
            bacterias_plot.append(b)
    color_map = [ plt.get_cmap("tab10")(i) for i in range(13)]
    color_map[10] = "peru"
    color_map[11] = "tab:orange"
    color_map[12] = "lawngreen"
    for i, b in enumerate(bacterias_plot):
        plots.append(
            plt.plot(
                stage_list,
                [ -log10(p_val[b]).real for p_val in p_val_dicts.values() ],
                marker = "o",
                markersize = 10,
                color = color_map[i]
            )[0]
        )
    # ax.set_yscale("log")
    # import matplotlib.ticker
    # ax.yaxis.set_major_locator(
    #     matplotlib.ticker.LogLocator( numticks=14 )
    # )
    # ax.yaxis.set_minor_locator(
    #     matplotlib.ticker.LogLocator(
    #         numticks=14,
    #         subs=( 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2 )
    #     ) 
    # )
    # ax.set_yticks( [10**i for i in range(0,13)] )
    # ax.set_yticklabels( [1]+[f"1e-{i}" for i in range(1,13)], minor=False )
    # ax.set_yticklabels( [], minor=True )

    ax.axhline(y=-log10(0.005).real, color=ColorCode.BLACK.value, lw=FIG_LINE_WIDTH)
    ax.set_ylabel("$-\log_{10}p$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        plots, bacterias_plot, 
        bbox_to_anchor=(1.02,0),
        loc = 'lower left',
        frameon = False
    )
    plt.savefig(f"for_paper_output/p_val_in3stage.pdf", bbox_inches='tight' )

def Output_averageAbundanceInSomeStages( stages, bacterias, stageToBactToAbdAvg ):
    import matplotlib.pyplot as plt
    from math import log10
    b_list = sorted( bacterias, reverse=True, key=lambda b:stageToBactToAbdAvg[stages[-1]][b] ) 
    plt.figure()
    plots = []
    for b in b_list:
        plots.append(
            plt.plot(
                stages,
                [ log10(stageToBactToAbdAvg[s][b]).real for s in stages ],
                marker = "o"
            )[0]
        )
    plt.legend(
        plots, b_list,
        bbox_to_anchor=(1.02,0),
        loc = 'lower left'
    )
    plt.savefig(f"for_paper_output/abundance_average_in{len(stages)}stages.png", bbox_inches='tight' )

def Output_Scatter_plot( 
        filename, 
        parm_x_label, parm_y_label, 
        bact_to_parm_x, bact_to_parm_y,
        bacteria_focused, bacterias,
        bact_to_abundance_mean,
        is_log_scale_y,
        log_ticks_reversed_y,
        draw_line=True
    ):
    import matplotlib.pyplot as plt
    from adjustText import adjust_text
    from utility import ColorCode
    colors = [ ColorCode.RED.value if b in bacteria_focused else ColorCode.GRAY.value for b in bacterias ]
    from math import log10

    min_abundance = min([bact_to_abundance_mean[b] for b in bacterias])
    sizes = [ AbundanceToCircleSize( bact_to_abundance_mean[b], min_abundance, 10, 100 ) for b in bacterias]

    xs = []
    ys = []
    for b in bacterias:
        x = bact_to_parm_x[b]
        y = bact_to_parm_y[b]
        xs.append(x)
        ys.append(y)
        # plt.annotate( b, (x,y) )
    
    _, ax = plt.subplots()
    ax.scatter( xs, ys, c = colors, s=sizes )

    # x=yの直線
    # ax.axline( [0,0], [7,7], color=ColorCode.BLACK )
    
    ax.set_xlabel( parm_x_label )
    ax.set_ylabel( parm_y_label )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    if is_log_scale_y:
        ax.set_yscale("log")
        if log_ticks_reversed_y:
            import matplotlib.ticker
            ax.yaxis.set_major_locator(
                matplotlib.ticker.LogLocator( numticks=7 )
            )
            ax.yaxis.set_minor_locator(
                matplotlib.ticker.LogLocator(
                    numticks=7,
                    subs=( 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2 )
                ) 
            )
            ax.set_yticks( [10**i for i in range(0,6)] )
            ax.set_yticklabels( [1]+[f"1e-{i}" for i in range(1,6)], minor=False )
            ax.set_yticklabels( [], minor=True )
    
            # p値の場合  
            if draw_line:
                ax.axhline(y=1/0.005, color=ColorCode.BLACK.value)
    elif draw_line:
        # p値の有意水準(logscaleでない場合)
        from math import log10
        ax.axhline( y=-log10(0.005).real, color=ColorCode.BLACK.value, lw=FIG_LINE_WIDTH )

    labels = []
    for i in range(len(bacterias)):
        x = xs[i]
        y = ys[i]
        b = bacterias[i]
        labels.append( ax.text(x,y,b) )

    adjust_text( labels , force_text=0.5, arrowprops = dict(arrowstyle="-", color=ColorCode.BLACK.value, alpha=0.5) )
    plt.savefig( filename )


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


def PPR_Diff( nx_g1, nx_g2 ):
    import networkx as nx
    ppr_diff_dict = dict()
    for focus_v in nx_g1.nodes():
        ppr_dict_1 = nx.pagerank( nx_g1, alpha=0.85, personalization={focus_v:1.0} )
        ppr_dict_2 = nx.pagerank( nx_g2, alpha=0.85, personalization={focus_v:1.0} )
        ppr_diff_dict[focus_v] = sum( abs(ppr_dict_2[v]-ppr_dict_1[v]) for v in nx_g1.nodes() )
    return ppr_diff_dict
    

def PPR_on_multipleGraph():
    pass

import numpy as np
def GraphLaplacian( V, edgeToWeight, convertBoolean=False ):
    deg = dict()
    for v1 in V:
        deg[v1] = 0
        for v2 in V:
            if (v1,v2) in edgeToWeight:
                deg[v1] += edgeToWeight[v1,v2]
            elif (v2,v1) in edgeToWeight:
                deg[v1] += edgeToWeight[v2,v1]
    array2d =[]
    for v1 in V:
        array = []
        for v2 in V:
            if v1==v2:
                array.append( deg[v1] )
            elif (v1, v2) in edgeToWeight:
                array.append( -edgeToWeight[v1,v2] )
            elif (v2, v1) in edgeToWeight:
                array.append( -edgeToWeight[v2,v1] )
            else:
                array.append( 0 )
        array2d.append( array )
    return np.array(array2d)
 
def NormalizedLaplacian( V, edgeToWeight ):
    degree = dict()
    for v1 in V:
        degree[v1] = 0
        for v2 in V:
            val = 0
            if (v1,v2) in edgeToWeight:
                val = edgeToWeight[v1,v2]
            elif (v2, v1) in edgeToWeight:
                val = edgeToWeight[v2,v1]
            degree[v1] += val


    nl_2darray = []
    for v1 in V:
        tmp_array = []
        for v2 in V:
            val = 0
            if v1==v2 and degree[v1] != 0:
                val = 1
            elif v1!=v2:
                if (v1,v2) in edgeToWeight:
                    val = - edgeToWeight[v1,v2] / math.sqrt( degree[v1] * degree[v2] )
                elif (v2, v1) in edgeToWeight:
                    val = -edgeToWeight[v2,v1] /  math.sqrt( degree[v1] * degree[v2] )
            tmp_array.append( val )
        nl_2darray.append(tmp_array)
    return np.array(nl_2darray)
                    

def SpectralClustering( V, edgeToWeight, cluster_num ):
    L = GraphLaplacian( V, edgeToWeight)
    
    eig_vals,eig_vecs = np.linalg.linalg.eig(L)

    eig_vecs = eig_vecs[:,np.argsort(eig_vals)]
    eig_vals = eig_vals[np.argsort(eig_vals)]

    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=cluster_num)
    kmeans.fit(eig_vecs[:,1:cluster_num])


def Calc_SDD( V, stageToEdgeToWeight, focusStage, spectreDim ):
    stages = list( stageToEdgeToWeight.keys() )
    stageToNodeToSpectre = dict()
    for s in stages:
        stageToNodeToSpectre[s] = CalcNodeToSpectre(V, stageToEdgeToWeight[s], spectreDim )

    stageToNodeToDistanceVec = dict()
    df_dict = dict()
    for s in stages:
        # if stageToK['Healthy']==26:
        #     print(f"======={s}===========")
        #     for v in V:
        #         print(f"{v}: ",end="")
        #         for val in nodeToSpectre[v]:
        #             print(f"\t{val:.3f}",end="")
        #         print()
        #     if s=='Stage\_0':
        #         exit()
        nodeToSpectre = stageToNodeToSpectre[s]
        nodeToDistanceVec = dict()
        for v1 in V:
            nodeToDistanceVec[v1] = [np.linalg.norm( nodeToSpectre[v2] - nodeToSpectre[v1] ) for v2 in V ]
        stageToNodeToDistanceVec[s] = nodeToDistanceVec
        df_dict[f"{v1} in {s}"] = nodeToDistanceVec

        '''
        print(f"========{s}============")
        for i, val in enumerate(eig_vals):
            print(f"{val}: {list(eig_vecs[:,i])}")
        print(f"10th eig_val: {eig_vals[9]}")
        '''
    
    import pandas as pd
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("display.precision", 6)

    # for s in stages:
        # print( pd.DataFrame(stageToNodeToDistanceVec[s]), file = open( f"for_theory_output/normalized/sdd/spectralDistance_{s}.txt", "w" ) )


    Y = focusStage
    if Y==stages[0]:
        X = stages[1]
    else:
        X = stages[0]    
    # vToMinusSum = dict()
    # vToPlusSum = dict()
        
    vToSDD = dict()
    for v in V:
        diff_vec = [ y-x for y,x in zip(stageToNodeToDistanceVec[Y][v], stageToNodeToDistanceVec[X][v]) ]
        # vToPlusSum[v] = sum( [x for x in diff_vec if x>0] )
        # vToMinusSum[v] = sum( [x for x in diff_vec if x<0] )
        vToSDD[v] = sum( [abs(x) for x in diff_vec] )

    return vToSDD

def Calc_SpectralDiff(V, stageToEdgeToWeight, focusStage, spectreDim ):
    stages = list( stageToEdgeToWeight.keys() )
    stageToNodeToSpectre = dict()
    for s in stages:
        stageToNodeToSpectre[s] = CalcNodeToSpectre(V, stageToEdgeToWeight[s], spectreDim )

    Y = focusStage
    if Y==stages[0]:
        X = stages[1]
    else:
        X = stages[0]

    vToSpectreDiff = dict()
    
    
    for v in V:
        # diff_vec = [ y-x for y,x in zip(stageToNodeToSpectre[Y][v], stageToNodeToSpectre[X][v]) ]
        # vToSpectreDiff[v] = sum( [abs(x) for x in diff_vec] )
        vToSpectreDiff[v] = min( np.linalg.norm( stageToNodeToSpectre[Y][v] - stageToNodeToSpectre[X][v] ), np.linalg.norm( stageToNodeToSpectre[Y][v] + stageToNodeToSpectre[X][v] ) )
    return vToSpectreDiff

def CalcNodeToSpectre(V, edgeToWeight, spectreDim, normalized=True ):
    if normalized:
        L = NormalizedLaplacian( V, edgeToWeight )
    else:
        L = GraphLaplacian( V, edgeToWeight, convertBoolean=False )

    eig_vals, eig_vecs = np.linalg.linalg.eig(L)
    eig_vecs = eig_vecs[:,np.argsort(eig_vals)]
    eig_vals = eig_vals[np.argsort(eig_vals)]

    # for v in eig_vecs:
    #     for val in v:
    #         print(val,end="\t")
    #     print()

    use_eig_vecs = [] 
    for i, val in enumerate(eig_vals):
        if i==0: #最小の固有値に対応する固有ベクトルはとばす
            continue
        use_eig_vecs.append(eig_vecs[:,i])
        if len(use_eig_vecs)>= spectreDim:
            break
    use_eig_vecs = np.array(use_eig_vecs)

    nodeToSpectre = dict()
    for i, v in enumerate(V):
        nodeToSpectre[v] = np.array(use_eig_vecs[:,i])
    return nodeToSpectre

    
# ネットワークが連結であること前提
def LeverageScore(V, edgeToWeight, spectreDim ):
    # nodeToSpectre = CalcNodeToSpectre( V, edgeToWeight, spectreDim, normalized=False )
    nodeToSpectre = CalcNodeToSpectre( V, edgeToWeight, spectreDim )
    nodeToLeverageScore = dict()
    for n in V:
        nodeToLeverageScore[n] = np.linalg.norm( nodeToSpectre[n] )
    return nodeToLeverageScore


def main():
    STEP3_SKIP = True
    CONTRACT_SKIP = False
    PICKUP_SKIP = False

    SKIP_TO_STEP7 = False
    args = sys.argv
    if len(args)<=3 :
        print("arg size must be 3", file=sys.stderr)
        return  

    for i in range(2):
        if not path.isfile(args[1]):
            v = ["table","sample"]
            print(f"{utility.ToOrdinal(i)} arg is not file path (it must be {v[i]} file path)", file=sys.stderr )
            return
    
    tableFilePath = args[1]
    sampleFilePath = args[2]
    focusStage = args[3]

    import pandas as pd
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("display.precision", 6)

    bacterias, stages,  bactToStageToCnts, stageToSamples = ReadCntTable_Whole( tableFilePath, sampleFilePath )
    # with open('sample_to_category_H_S0.tsv', 'w') as f:
    #     f.write('Sample\tCategory\n')
    #     for s in ['Healthy', 'Stage_0']:
    #         for sample in stageToSamples[s]:
    #             f.write(f"{sample}\t{s}\n")
    # with open("genus_count_table_H_S0.tsv", "w") as f:
    #     for s in ['Healthy', 'Stage_0']:
    #         for i, sample in enumerate(stageToSamples[s]):
    #             if not (s=='Healthy' and i==0) :
    #                 f.write('\t')    
    #             f.write(sample)
    #     f.write('\n')
    #     for b in bacterias:
    #         f.write(b)
    #         for s in ['Healthy', 'Stage_0']:
    #             for cnt in bactToStageToCnts[b][s]:
    #                 f.write(f'\t{cnt}')
    #         f.write('\n')

    #bacterias, stages = ReadCntTable( tableFilePath )
    
    bactToStageToAbds = Step1( bacterias, stages, bactToStageToCnts )
    print("step1 finished")
    #大腸がんデータに対する特殊処理 ツール化する際に消す
    stages = ["Healthy", focusStage]
    # stages = ["Healthy", "Multiple_polyps","Stage_0","Stage_I_II","Stage_III_IV"]

    if args[3] not in stages:
        print(f"Third arg:{args[3]} is not found from stage column in sampleFile.", file=sys.stderr)
        return

    if stages[0]==focusStage:
        baseStage = stages[1]
    else:
        baseStage = stages[0]

    _, before_contract_p_val = Step2( MWU_P_TH, MWU_RANK_TH, baseStage, focusStage, bactToStageToAbds )
    #print("step2 finished")
    
    if STEP3_SKIP:
        print("step3 skipped")
    else:
        Step3( stages, bacterias, bactToStageToCnts, stageToSamples )
        print("step3 finished")
    
    stageToEdgeToWeight, bacteria_effective, stageToVTupleToCorr = Step4( EDGE_THRESHOLD, stages, bacterias)

    print("step4 finished")

    # Output_GraphInfo_forTeX(
    #     "for_paper_output/graph_info/graph_info_forTeX_allStage.txt",
    #     bacteria_effective,
    #     stageToEdgeToWeight
    # )
    # exit()

    abundance_mean = CalcAllAbundanceMean( bactToStageToAbds )

    sup_category = ReadCategory( "category.tsv" )

    '''平均abundanceのフィルタ
    temp_vec = []
    for b in bacteria_effective:
        if b not in sup_category.keys():
            print(f"{b} removed because not in category.tsv")
        elif abundance_mean[b] < 1e-9:
            print(f"{b} removed because abundance < 1e-9")
        else:
            temp_vec.append(b)
    bacteria_effective = temp_vec
    '''
    
    # if CONTRACT_SKIP:
    #     bacteria_contracted = bacteria_effective
    #     stageToEdgeToWeight_contracted = stageToEdgeToWeight
    #     print("step5 skipped")
    # else:
    stageToBactToCluster = Step5_1( bacteria_effective, stageToEdgeToWeight )

    parentToContractSet = Step5_2( bacteria_effective, sup_category, abundance_mean, stageToBactToCluster )

    bacteria_contracted, stageToEdgeToWeight_contracted = Step5_3( parentToContractSet, stageToEdgeToWeight )
    print("step5 finished")

    bactToStageToAbds = { b:bactToStageToAbds[b] for b in bacteria_contracted }
    bacteria_focused, p_val = Step2( MWU_P_TH, MWU_RANK_TH, baseStage, focusStage, bactToStageToAbds )
    print("step2 finished")

    stageTo_VE_tuple = Step6( PICKUP_SKIP, bacteria_focused, bacteria_contracted, stageToEdgeToWeight_contracted )
    if PICKUP_SKIP:
        print("step6 skipped")
    else :
        print("step6 finished")

    V = stageTo_VE_tuple[stages[0]][0]
    stageToEdgeToWeight_picked = {s:stageTo_VE_tuple[s][1] for s in stageTo_VE_tuple.keys()}
    '''辺重みを0.4以下も含める処理
    _, stageToEdgeToWeight_withUnderTH = Step5_3( parentToContractSet, stageToVTupleToCorr )
    for s in stages:
        for v1 in V:
            for v2 in V:
                if v1==v2:
                    continue
                e = (v1,v2)
                if e in stageToEdgeToWeight_withUnderTH[s]:
                    if stageToEdgeToWeight_withUnderTH[s][e] > 0.0:
                        stageTo_VE_tuple[s][1][e] = stageToEdgeToWeight_withUnderTH[s][e]
                else:
                    if stageToEdgeToWeight_withUnderTH[s][e[1],e[0]] > 0.0:
                        stageTo_VE_tuple[s][1][e] = stageToEdgeToWeight_withUnderTH[s][e[1],e[0]]
    '''
    '''辺のないところに0.1の辺を追加
    for s in stages:
        for v1 in V:
            for v2 in V:
                if v1==v2:
                    continue
                e = (v1,v2)
                re = (v2,v1)
                if e not in stageToEdgeToWeight_picked[s] and re not in stageToEdgeToWeight_picked[s]:
                    stageToEdgeToWeight_picked[s][e] = 0.01
    '''

    #     df = pd.DataFrame( data= stageTo_VE_tuple )
    #     df.to_csv('mid_datas/step6result.tsv', sep="\t")

    #     df = pd.DataFrame( data = stageToEdgeToWeight_withUnderTH )
    #     df.to_csv('mid_datas/stageToEdgeToWeight_withUnder04.tsv', sep="\t")

    #     with open ('mid_datas/abundance_mean.tsv', 'w') as f:
    #         for v in abundance_mean.keys():
    #             f.write(f"{v}\t{abundance_mean[v]}\n")
    
    # if SKIP_TO_STEP7:
    #     stageTo_VE_tuple = dict()
    #     f = open('mid_datas/step6result.tsv',"r")
    #     for line in f:
    #         words = line.split()
    #         if len(words) < 4:
    #             continue
    #         words

    #     abundance_mean = dict()
    #     abundance_mean_f =  open('mid_datas/abundance_mean.tsv', 'r')
    #     for line in abundance_mean_f:
    #         tuple = line.split()
    #         abundance_mean[tuple[0]] = float(tuple[1])
    #     stageToEdgeToWeight_withUnderTH = {stages[0]:dict(), stages[1]:dict()}
    #     f = open('mid_datas/stageToEdgeToWeight_withUnder04.tsv')
    #     for line in f:
    #         words = line.split()
    #         if len(words) < 4:
    #             continue
    #         edge = (words[0], words[1])
    #         stageToEdgeToWeight_withUnderTH[stages[0]][edge] =  words[2]
    #         stageToEdgeToWeight_withUnderTH[stages[1]][edge] =  words[3]

    #     V = stageTo_VE_tuple[stages[0]][0]
        

    # _, p_val_in_S12 = Step2( MWU_P_TH, MWU_RANK_TH, stages, "Stage_I_II", bactToStageToAbds )
    # _, p_val_in_S34 = Step2( MWU_P_TH, MWU_RANK_TH, stages, "Stage_III_IV", bactToStageToAbds )
    
    abundance_mean_inGraph = {b:abundance_mean[b] for b in V}

    stageToBactToAbdAvg = { s:{ b: sum(bactToStageToAbds[b][s]) / len(bactToStageToAbds[b][s]) for b in V } for s in stages}
    sym_diff, nesh, weighted_nesh, stageToBactToDegree = Step7( MWU_RANK_TH, CONTRACT_SKIP, PICKUP_SKIP, 
        stageTo_VE_tuple, 
        focusStage, 
        p_val, 
        bacteria_focused, 
        abundance_mean_inGraph, 
        stageToBactToAbdAvg
    )
    print("step7 finished")

    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = 'Bahnschrift'
    plt.rcParams["font.size"] = 10
    # "Arial" "Bahnschrift" 'Times New Roman' 'Yu Gothic'
    
    
    '''ネットワークをテキストで出力
    displayed = {v:False for v in V}
    for s in stageToEdgeToWeight_picked.keys():
         with open(f"for_houkoku_output/network_{s}.tsv", "w")as f:
            for e in stageToEdgeToWeight_picked[s].keys():
                f.write(f"{e[0]}\t{e[1]}\t{stageToEdgeToWeight_picked[s][e]}\n")
                displayed[e[0]] = True
                displayed[e[1]] = True
            for v in V:
                if not displayed[v]:
                    f.write(f"{v}\n")
    # 隣接行列ver
    for s in stageToEdgeToWeight_picked.keys():
         with open(f"for_houkoku_output/adjacency_matrix_{s}.tsv", "w")as f:
            for v in V:
                f.write(f"\t{v}")
            f.write("\n")
            for v1 in V:
                f.write(f"{v1}")
                for v2 in V:
                    if v1==v2:
                        f.write("\t1")
                    elif (v1,v2) in stageToEdgeToWeight_picked[s]:
                        w = stageToEdgeToWeight_picked[s][v1,v2]
                        f.write(f"\t{w}")
                    elif (v2,v1) in stageToEdgeToWeight_picked[s]:
                        w = stageToEdgeToWeight_picked[s][v2,v1]
                        f.write(f"\t{w}")
                    else:
                        f.write("\t0")
                f.write("\n")
    '''
    
    '''abundanceをテキスト出力
    for s in stages:
        with open(f"for_houkoku_output/abundance_{s}_H-{focusStage}.txt", "w") as f:
            for sample in stageToSamples[s]:
                f.write(f"\t{sample}")
            f.write("\n")
            for v in V:
                f.write(v)
                for abd in bactToStageToAbds[v][s]:
                    f.write(f"\t{abd}")
                f.write("\n")

    '''

    '''#p値をテキスト出力
    with open("for_houkoku_output/p-val.txt","w") as f:
        f.write("細菌名(属)\tp値\n")
        for label, val in value_sorted_items:
            f.write(f"{label}\t{val}\n")
    '''
    
    '''#p値の棒グラフ出力 
    from math import log10
    minuslog10pValDic = { k:-log10(before_contract_p_val[k]).real for k in before_contract_p_val.keys() if before_contract_p_val[k] < 0.01 }
    value_sorted_items = sorted(list(minuslog10pValDic.items()), key=lambda x:x[1] )
    labels, values = utility.UnpackListOfTuple( value_sorted_items )

    from utility import ColorCode
    colors = [ ColorCode.RED.value if l in bacteria_focused else ColorCode.GRAY.value for l in labels ]

    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.barh( labels, values, color=colors)
    ax.axvline( x = -log10(0.005).real, color=ColorCode.BLACK.value, lw = FIG_LINE_WIDTH )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("$-\log_{10}p$")
    plt.savefig("for_paper_output/p-val_bar_linear.pdf", bbox_inches='tight', pad_inches=0 )
    '''


    '''#log目盛バージョン
    minuspValDic = { k:1/p_val[k] for k in p_val.keys() if p_val[k] < 0.01 }
    value_sorted_items = sorted(list(minuspValDic.items()), key=lambda x:x[1] )
    labels, values = utility.UnpackListOfTuple( value_sorted_items )from utility import ColorCode
    colors = [ ColorCode.RED.value if l in bacteria_focused else ColorCode.GRAY.value for l in labels ]
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(
        matplotlib.ticker.LogLocator( numticks=5 )
    )
    ax.xaxis.set_minor_locator(
        matplotlib.ticker.LogLocator(
            numticks=5,
            subs=( 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2 )
        )
    )
    ax.set_xticks( [10**i for i in range(2,6)] )
    ax.set_xticklabels( [f"1e-{i}" for i in range(2,6)], minor=False )
    ax.set_xticklabels( [], minor=True )

    ax.set_yticks( labels, fontsize=5 )
    
    ax.axvline(x=1/0.005, color=ColorCode.BLACK.value)
    plt.savefig("for_paper_output/p-val_barh.png",bbox_inches='tight', pad_inches=0)
    '''
    
    #Output_ContractInfo( parentToContractSet, sup_category )
   
    '''#p値変化を折れ線で出力
    research_stages = ["Stage_0","Stage_I_II","Stage_III_IV"]
    p_val_dicts = dict()
    for s in research_stages:
        _, p_val_dicts[s] = Step2( MWU_P_TH, MWU_RANK_TH, baseStage, s, bactToStageToAbds )
    Output_pValInSomeStages( V, p_val_dicts, MWU_P_TH )
    '''
    
    '''#average abundance変化を折れ線で出力
    research_stages = ["Healthy","Stage_0","Stage_I_II", "Stage_III_IV"]
    Output_averageAbundanceInSomeStages( research_stages, V, stageToBactToAbdAvg)
    exit()
    '''
    
    '''縮約前ネットワークの基本情報出力
    Output_GraphInfo_forTeX(
        "for_paper_output/graph_info/graph_info_forTeX.txt",
        bacteria_effective,
        stageToEdgeToWeight
    )
    Output_GraphInfo(
        "for_paper_output/graph_info/graph_info.txt",
        bacteria_effective,
        stageToEdgeToWeight
    )
    '''

    '''#縮約後のネットワークの基本情報出力
    Output_GraphInfo_forTeX(
        f"for_paper_output/graph_info/graph_contracted_info_{focusStage}_th=0.txt",
        bacteria_contracted,
        stageToEdgeToWeight_contracted
    )
    '''
    '''最終的なネットワークの基本情報出力
    Output_GraphInfo_forTeX(
        f"for_paper_output/graph_info/graph_picked_info_{focusStage}.txt",
        V,
        stageToEdgeToWeight_picked
    )
    '''
    
    '''#縮約情報の表出力
    minusSize_parent_ContractSet_tupleV = [ (-len(contractSet),parent,contractSet) for parent,contractSet in parentToContractSet.items() ]
    with open("for_houkoku_output/contract_info_1214_th=0.txt","w")  as f:
        f.write("細菌群の & 類似細菌群(先頭{\\bf 太字}が代表菌) \\\\大きさ & \\\\ \\hline\\hline \n")
        single_bacts=[]
        for minusSize, parent, contractSet in sorted(minusSize_parent_ContractSet_tupleV):
            if -minusSize==1:
                single_bacts.append(parent)
                continue                
            f.write(f"{-minusSize} & {{\\bf {parent} }}")
            for b in sorted(contractSet):
                if b == parent:
                    continue
                f.write(f", {b}")
            f.write("\\\\ \\hline \n")
        f.write("==============一元化を経ずにbacteria_contractedに残っているもの==============\n")
        for b in sorted(single_bacts):
            f.write(f"{b},")
    '''
    
    '''指定した1菌の周辺についてグラフを出力
    ResearchOneBact_beforeAndAfterContract( 
        research_bact_name = "Obesumbacterium",
        whole_bacterias = bacteria_effective,
        parentToContractSet = parentToContractSet,
        stageToEdgeToWeight_before = stageToAdjacency_matrix,
        stageToEdgeToWeight_after = stageToEdgeToWeight_contracted,
        stageToBactToCluster = stageToBactToCluster
    )'''
    
    '''グラフ出力
    # 縮約前
    # stage_to_before_contract_VE_tuple = {
    #     s:(
    #         bacteria_effective, 
    #         { edge:weight for edge,weight in stageToEdgeToWeight[s].items() if edge[0] in bacteria_effective and edge[1] in bacteria_effective }
    #     )
    #     for s in stages
    # }
    # Output_NxGraph( "_before_contract",stage_to_before_contract_VE_tuple, focusStage, bacteria_focused, abundance_mean, draw_label=False, node_size_change=False, edge_width_change=False, sort_as_degree=False )
    
    # 縮約後
    stage_to_before_contract_VE_tuple = {
        s:(
            bacteria_contracted, 
            stageToEdgeToWeight_contracted[s]
        )
        for s in stages
    }
    Output_NxGraph( "_contracted",stage_to_before_contract_VE_tuple, focusStage, bacteria_focused, abundance_mean, draw_label=False, node_size_change=False, edge_width_change=False, sort_as_degree=False )
    
    # Output_ForCytoscape( G, bacteria_focused )
    # Output_ForNetShiftApp(  )
    
    # pickup後（最終）
    '''
    Output_NxGraph( "", stageTo_VE_tuple, sym_diff, bacteria_focused, abundance_mean_inGraph )

    ''' 差分ネットワーク
    degreeInFocusStage = dict()
    for v in V:
        degreeInFocusStage[v] = 0.0
    for e in stageToEdgeToWeight_picked[focusStage]:
        degreeInFocusStage[e[0]] += stageToEdgeToWeight_picked[focusStage][e]
        degreeInFocusStage[e[1]] += stageToEdgeToWeight_picked[focusStage][e]
        
    sorted_V = sorted( V, key=lambda v:sym_diff[v], reverse=True )
    node_size_dict = Convert_Feature_to_Size( sym_diff, [100, 1200] )
    Output_DifferGraph(
        sorted_V, 
        bacteria_focused, 
        stageToEdgeToWeight_picked,
        [node_size_dict[v] for v in sorted_V]
    )
    '''
    
    #SpectralDiistanceDiff
    ''' # plus_sum, minus_sum, all_sum, abs_sumそれぞれについて
    
    sort_key_features = [
        "abs_sum",
        #"all_sum"
    ]
    df_dict = {
        "minus_sum":minus_sum,
        "plus_sum":plus_sum,
        "all_sum":{ v:minus_sum[v]+plus_sum[v] for v in V},
        "abs_sum":{ v:plus_sum[v]-minus_sum[v] for v in V}
    }
    output_df = pd.DataFrame(df_dict)

    for val_name in sort_key_features:
        print(
            output_df.sort_values(val_name, ascending = False), 
            file=open(f"for_paper_output/spectralDistanceDiff_sortBy{val_name}.txt","w") 
        )

    '''
    
    '''# kをそれぞれの群で変えながら
    df_dict = dict()
    k_list_H = [13,24]
    k_list_S0 = [6,13,22]
    for k_H in k_list_H:
        for k_S0 in k_list_S0:
             = Calc_SpectralDistanceDiff( V, stageToEdgeToWeight, focusStage, {stages[0]:k_H,stages[1]:k_S0} )
            df_dict[k_H,k_S0] = { v:plus_sum[v]-minus_sum[v] for v in V }
    
    output_df = pd.DataFrame(df_dict)
    print(
        output_df.sort_values((24,22), ascending = False), 
        file=open(f"for_paper_output/spectralDistanceDiff_sortBy24_22.txt","w") 
    )
    '''
    
    # kToVToSDD = dict()
    # range_v = list(range(1,len(V),3))
    # range_v.append(25)
    # range_v = sorted(range_v)
    '''SDD
    sdd_df_dict = dict()
    for k_common in range(1,len(V)):
        vToSDD = Calc_SDD( V, stageToEdgeToWeight_picked, focusStage, k_common )
        sdd_df_dict[k_common] = vToSDD
    
    output_df = pd.DataFrame(sdd_df_dict)
    print(
        output_df.sort_values(1, ascending = False), 
        file=open(f"for_theory_output/normalized/sdd/sdd_same_k_sortByk={len(V)-1}_withunder04edges.txt","w") 
    )
    '''
    '''tex用 
    print("\\hline 頂点 & k=", end="")
    for k in range_v:
        print(k, end="")
        if k<=len(V)-1:
            print(" & ", end="")
    print("\\\\ \\hline")
    for v in sorted(V, key=lambda x:kToVToSDD[range_v[-2]][x], reverse=True):
        print(f"{v}", end="")
        for k in range_v:
            print(f" & ${kToVToSDD[k][v]:.3f}$", end="")
        print("\\\\ \\hline")
    '''
    

    '''# SpectreDiff
    spectreDiff_df_dict = dict()
    spectre_dim_list = list(range(1, len(V)))
    for spectre_dim in spectre_dim_list:
        vToSpectreDiff = Calc_SpectralDiff( V, stageToEdgeToWeight_picked, focusStage, spectre_dim)
        spectreDiff_df_dict[spectre_dim] = vToSpectreDiff
        
    output_df = pd.DataFrame(spectreDiff_df_dict)
    print(
        output_df.sort_values(spectre_dim_list[0], ascending = False), 
        file=open(f"for_theory_output/normalized/spectreDiff/spectralDiff_sortByDim={spectre_dim_list[0]}.txt","w") 
    )
    '''
    

    '''結果のtex用出力
    for_tex_df_dict = {
        "QNetDiff": sym_diff,
        "NetShift": nesh
        # "p-value(abundance_elevated)":{b:p_val[b] for b in bacterias},
        # "average_abundance":{ b:abundance_mean[b] for b in bacterias},
    }
    For_tex_output( "for_paper_output/QNetDiff_Netshift_result_allgray_fortex.txt", for_tex_df_dict, V, {})
    '''

    # with open("for_slide_output/QNet.csv","w") as of:
    #     for v in sorted(V,key=lambda x:sym_diff[x], reverse=True):
    #         print(f" {v},{sym_diff[v]:.3f} ", file=of)
    

    
    import networkx as nx
    '''PPR
    nx_g1 = GetNx_g(V, stageToEdgeToWeight_picked[stages[0]])
    nx_g2 = GetNx_g(V, stageToEdgeToWeight_picked[stages[1]])
    dicts = dict()

    ppr_diff_dict = dict()
    personalize_rate = 1.0
    for focus_v in V:
        # print(f"========== focus_v: {focus_v} ==============")
        # print("in Healthy", end="" )
        ppr_dict1 = nx.pagerank( nx_g1, alpha=0.85, personalization={ v: (1-personalize_rate)/len(V) + (personalize_rate if v==focus_v else 0) for v in V } )
        # for val in ppr_dict.values():
        #     print(f"\t{val:.3f}", end="")
        dicts[f"{focus_v};H"] = ppr_dict1

        # print("\nin Stage_0", end="" )
        ppr_dict2 = nx.pagerank( nx_g2, alpha=0.85, personalization={ v: (1-personalize_rate)/len(V) + (personalize_rate if v==focus_v else 0) for v in V } )
        # for val in ppr_dict.values():
        #     print(f"\t{val:.3f}", end="")
        # print()
        dicts[f"{focus_v};S0"] = ppr_dict2

        dicts[f"{focus_v};S0-H"] = {v: ppr_dict2[v]-ppr_dict1[v] for v in V}
        ppr_diff_dict[focus_v] = sum( abs(ppr_dict2[v]-ppr_dict1[v]) for v in V)


    # output_df = pd.DataFrame( data=dicts )
    # print(
    #     #output_df.sort_values(f"{stages[1]} - {stages[0]}", ascending = False), 
    #     output_df,
    #     file=open(f"for_houkoku_output/ppr/all_ppr_witheps=001.txt","w") 
    # )
    # output_df = pd.DataFrame({"ppr_diff":ppr_diff_dict})
    # print(
    #     output_df.sort_values("ppr_diff",ascending = False), 
    #     file=open(f"for_houkoku_output/ppr/ppr_diff_witheps=001.txt","w") 
    # )

    '''
    
    
    '''page rank
    stage_to_pagerank_dict = dict()
    for s in stages:
        nx_g = GetNx_g(V, stageToEdgeToWeight_picked[s])
        stage_to_pagerank_dict[s] = nx.pagerank( nx_g, alpha=0.85 )
    # output_df = pd.DataFrame(
    #     {
    #         f"{stages[1]} - {stages[0]}":{b:abs(stage_to_pagerank_dict[stages[1]][b]-stage_to_pagerank_dict[stages[0]][b]) for b in V},
    #         f"pagerank in {stages[0]}":stage_to_pagerank_dict[stages[0]], 
    #         f"pagerank in {stages[1]}":stage_to_pagerank_dict[stages[1]],
    #         #f"abs diff ":{b:abs(stage_to_pagerank_dict[stages[1]][b]-stage_to_pagerank_dict[stages[0]][b]) for b in V},
    #     } 
    # )
    # print(
    #     output_df.sort_values(f"{stages[1]} - {stages[0]}", ascending = False), 
    #     file=open(f"for_houkoku_output/pagerank.txt","w") 
    # )
    '''

    '''# 散布図
    x_param_label_dict = {
        "sym_diff":sym_diff,
    #    "weighted NESH score":weighted_nesh,
    #    "degree in Stage_0":stageToBactToDegree[stages[1]]
    #    "NetShift score":nesh
    }
    from math import log10
    y_param_label_dict = {
    #    "$-\log_{10}p$":{ b: -log10(p_val[b]).real for b in V },
    #    "SDD(k=25)":kToVToSDD[25]
        "NetShift":nesh
    #    "pagerank diff":{b:abs(stage_to_pagerank_dict[stages[1]][b]-stage_to_pagerank_dict[stages[0]][b]) for b in V}
    }
    
    for x_param_label, x_param in x_param_label_dict.items():
        for y_param_label, y_param in y_param_label_dict.items():
            # bacterias = sorted(bacterias, key=lambda b:param_x[b], reverse=True )[:30] #その指標の上位30菌のみ
            y_label = y_param_label
            if y_param_label == "$-\log_{10}p$":
                y_label = "-log10p"
            Output_Scatter_plot( 
                f"for_paper_output/{x_param_label}-{y_label}_scatter.svg", 
                x_param_label, y_param_label, 
                x_param, y_param,
                bacteria_focused, V,
                abundance_mean,
                is_log_scale_y=False,
                log_ticks_reversed_y=False,
                draw_line=False       
            )
    '''
    
    '''グラフの特徴量の棒グラフ
    # label_to_param_v = dict()
    # index = 0
    # bacterias = sorted(bacterias, key=lambda b:degree_diff[b], reverse=True )[:30] #degree_diffの上位30菌のみ

    # for label in param_label_v:
    #     label_to_param_v[label] = [ param_dict_v[index][b] for b in bacterias ]
    #     index += 1
    # Output_MultiBar( 
    #     filepath = for_paper_path + "multibar.png", 
    #     keys = bacterias,
    #     param_labels = param_label_v,
    #     label_to_param_v = label_to_param_v,
    #     xtick_list = [ '*'+b if contract_remain_bacterias!=None and b not in contract_remain_bacterias else b  for b in bacterias]
    # )'''
    
    '''NetMoss用 abundanceファイル出力
    for s in stages:
        with open(f"for_netmoss_output/abundance_{s}.txt","w") as f:
            for i in range(len(stageToSamples[s])):
                f.write(f"\tsample_{i}")
            f.write("\n")
            for v in V:
                f.write(f"{v}\t")
                for cnt in bactToStageToCnts[v][s]:
                    f.write(f"{cnt}\t")
                f.write("\n")
    '''
    
    

    '''連結成分
    import networkx as nx
    nx_g1 = GetNx_g( V, stageToEdgeToWeight_picked[stages[0]] )
    nx_g2 = GetNx_g( V, stageToEdgeToWeight_picked[stages[1]] )
    for c in nx.connected_components(nx_g1):
        print(f"{len(c)} & ",end="")
        for b in c:
            print(b, end=", " if b != list(c)[-1] else "")
        print("\\\\ \\hline")
    '''

    '''レバレッジスコア
    kToNodeToLeverageDiff = dict()
    for spectreDim in range(1,len(V)):
        stageToNodeToLeverageScore = dict()
        for s in stages:
            stageToNodeToLeverageScore[s] = LeverageScore(V, stageToEdgeToWeight_picked[s], spectreDim )
        kToNodeToLeverageDiff[spectreDim] = {n:abs(stageToNodeToLeverageScore[stages[1]][n]-stageToNodeToLeverageScore[stages[0]][n]) for n in V }
        
        df = pd.DataFrame( 
            {
                f"Leverage score in {stages[0]}": stageToNodeToLeverageScore[stages[0]],
                f"Leverage score in {stages[1]}": stageToNodeToLeverageScore[stages[1]],
                "Leverage score diff": {n:abs(stageToNodeToLeverageScore[stages[1]][n]-stageToNodeToLeverageScore[stages[0]][n]) for n in V } 
            } 
        )
        print( df.sort_values("Leverage score diff", ascending=False), file=open(f"for_theory_output/normalized/Leverage/with04_LeverageScore_dim={spectreDim}","w"))
    df = pd.DataFrame( kToNodeToLeverageDiff )
    print( df.sort_values( 1, ascending=False ), file=open(f"for_theory_output/normalized/Leverage/with04_LeverageDiffs.txt","w") )
    '''
    


    '''# 各指標の表出力
    '''

    df_dict = {
    #    f"degree_diff_({focus_stage}-{non_focus_stage})": degree_diff,
    #    "edge_weight_jaccard_coefficient": jaccard,
    #    "weighted_nesh_score":weighted_nesh,
        "QNetDiff": sym_diff,
    #    "pagerank_diff": {b:abs(stage_to_pagerank_dict[stages[1]][b]-stage_to_pagerank_dict[stages[0]][b]) for b in V},
    #    "ppr_diff":ppr_diff_dict,
    #    f"pagerank in {stages[0]}": stage_to_pagerank_dict[stages[0]],
    #    f"pagerank in {stages[1]}": stage_to_pagerank_dict[stages[1]]
        "NetShift":nesh,
        # "p-value":{b:p_val[b] for b in V},
        # #f"isfocused(p-value_top{MWU_RANK_TH})":{b:b in bacteria_focused for b in bacterias},
        # #"average_abundance":{ b:abundance_mean[b] for b in V},
        # f"degree_in_{stages[0]}":stageToBactToDegree[stages[0]],
        # f"degree_in_{stages[1]}":stageToBactToDegree[stages[1]],
        # f"degree diff":{b:abs(stageToBactToDegree[stages[1]][b]-stageToBactToDegree[stages[0]][b]) for b in V},
        # f"abundance_average_in_{stages[0]}":stageToBactToAbdAvg[stages[0]],
        # f"abundance_average_in_{stages[1]}":stageToBactToAbdAvg[stages[1]]
        # #"SDD(k=25)":kToVToSDD[25],
    #    f"unique_edge_degree_in_{s0}":unique_degree[s0],
    #    f"unique_edge_degree_in_{s1}":unique_degree[s1]
    }
    output_df = pd.DataFrame(df_dict)
    
    output_folda_path = f"tool_output/"
    #output_folda_path = "tool_output/withUnder04edgeResult_"
    #output_folda_path = "tool_output/with_01edges"
    print( output_df.sort_values("QNetDiff", ascending=False) , file=open(output_folda_path+f"{stages[1]}-{stages[0]}_result_table.txt","w") )
    for val_name in df_dict.keys():
        if val_name == f"isfocused(p-value_top{MWU_RANK_TH})": 
            continue
        if val_name != "QNetDiff_score":
            continue
        print(output_df.sort_values(val_name, ascending = True if val_name == "p-value(abundance_elevated)" else False), file=open(output_folda_path+f"sort_by_{val_name}.txt","w") )

if __name__ == "__main__":
    main()