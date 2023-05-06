from cmath import e
from errno import EDEADLK
from inspect import CO_ASYNC_GENERATOR
from turtle import width
import networkx as nx
import matplotlib.pyplot as plt
from utility import ColorCode, similar_color_dict

def Show_net( nx_g, filename, node_to_size,  nodeGroup_to_color, edgecolor_v, node_shape='o', draw_label=True, edge_width_change=True, edge_weight_min=0.4, edge_weight_max=1.0 ):
    plt.figure( figsize = (10.0, 10.0) )
    groups_set = set( nx.get_node_attributes( nx_g, 'group').values() )
    lim = 1.85
    plt.xlim([-lim, lim])
    plt.ylim([-lim, lim])
    nodecolor_v = [ 
        nodeGroup_to_color[nx_g.nodes[n]['group']] 
        for n in nx_g.nodes
    ]

    # min_width = 1
    # max_width = 7
    min_width = 0.1
    max_width = 7
    
    width_grad = (max_width-min_width)/(edge_weight_max-edge_weight_min)
    # min_alpha = 0.3
    # max_alpha = 1.0
    min_alpha = 0.05
    max_alpha = 1.4
    
    alpha_glad = (max_alpha-min_alpha)/(edge_weight_max-edge_weight_min)

    if edge_width_change:
        edgewidth_v = [ min_width + (parm["weight"]-edge_weight_min)*width_grad for u, v, parm in nx_g.edges(data=True)]
        edgealpha_v = [ min_alpha + (parm["weight"]-edge_weight_min)*alpha_glad for u, v, parm in nx_g.edges(data=True)]
    else:
        edgewidth_v = [ 1 for u,v,_ in nx_g.edges(data=True)]
        edgealpha_v = [ 0.5 for u,v,_ in nx_g.edges(data=True)]
    #edgecolor_v = [ color_dict[code][1] for code in edge_colorCode_v ]

    # pos = nx.circular_layout(nx_g)
    pos = Nx_Circular_Pos( nx_g.nodes, 1.0 )
    nx.draw_networkx_nodes( 
        nx_g, pos, 
        node_shape = node_shape,
        node_color = nodecolor_v, 
    #    linewidths = 1.0, 
    #    edgecolors = nodeframecolor_v, 
        node_size=node_to_size
    )
    nx.draw_networkx_edges( 
        nx_g, pos, 
        alpha=edgealpha_v, 
        edge_color = edgecolor_v, 
        width = edgewidth_v 
    )
    
    if draw_label:
        nodes = nx_g.nodes
        node_num = len(nodes)
        left_half = []
        right_half = []
        for i,n in enumerate(nodes):
            if node_num/4 <= i and i < node_num/4*3:
                left_half.append(n)
            else:
                right_half.append(n)
        left_G = nx.Graph()
        left_G.add_nodes_from(left_half)
        right_G = nx.Graph()
        right_G.add_nodes_from(right_half)

        angle_per_node = 2*pi/node_num
        r = 1.15

        left_pos = dict()
        right_pos = dict()
        for i,n in enumerate(nodes):
            pos = ( r*cos( angle_per_node * i ), r*sin( angle_per_node * i ) )
            if len(nodes)/4 <= i and i < len(nodes)/4*3:
                left_pos[n] = pos
            else:
                right_pos[n] = pos 
        
        # MovePos( left_pos, 'Solobacterium', y = -0.030 )
        # MovePos( right_pos, 'Haemophilus', y = 0.025 )
        # MovePos( right_pos, 'Ruminococcus', y= -0.01 )
        # MovePos( left_pos, 'Veillonella', y=0.01 )

        nx.draw_networkx_labels( left_G, left_pos, horizontalalignment = 'right')
        nx.draw_networkx_labels( right_G, right_pos, horizontalalignment = 'left')
    # pos = Label_Circular_Pos( nx_g.nodes, 1.07 )
    # label_describe = nx.draw_networkx_labels( nx_g, pos, font_size=10.5, horizontalalignment = 'left')
    # for i, t in enumerate(label_describe.values()):
    #     angle = i*360/len(label_describe)
    #     if 90 < angle  <270:
    #         angle += 180
    #     t.set_rotation(angle)
    plt.savefig( filename, bbox_inches='tight', pad_inches=0 )

from math import cos, sin, pi
def Nx_Circular_Pos( nodes, r ):
    node_num = len(nodes)
    angle_per_node = 2*pi/node_num
    pos = dict()
    for i, n in enumerate(nodes) :
        pos[n] = ( r*cos( angle_per_node * i ), r*sin( angle_per_node * i ) ) 
        
    return pos


def Label_Circular_Pos(nodes,r):
    node_num = len(nodes)
    angle_per_node = 2*pi/node_num
    ONE_CHAR_WIDTH=0.02
    pos = dict()
    for i, n in enumerate(nodes) :
        if len(nodes)/4 <= i and i < len(nodes)/4*3:
            x_margin = len(n)*ONE_CHAR_WIDTH
        else:
            x_margin=0
        pos[n] = ( r*cos( angle_per_node * i ) - x_margin, r*sin( angle_per_node * i ) ) 
        
    return pos

def MovePos( pos_dict, key, x=0, y=0):
    pos_dict[key] = (pos_dict[key][0]+x, pos_dict[key][1]+y)