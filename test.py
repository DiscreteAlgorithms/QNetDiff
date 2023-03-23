# import matplotlib.font_manager as fm
# import pprint

# font_list = [f.name for f in fm.fontManager.ttflist]

# pprint.pprint(font_list)
# exit()

import numpy

def PrintMatrix( A ):
    for l in A:
        for i,x in enumerate(l):
            print(f"{x}",end="")
            if i != len(l)-1:
                print(" & ",end="")
        print("\\\\")
def GraphLaplacian( V, edgeToWeight ):
    deg = dict()
    for v1 in V:
        deg[v1] = 0
        for v2 in V:
            if (v1,v2) in edgeToWeight:
                deg[v1] += edgeToWeight[v1,v2]
            elif (v2,v1) in edgeToWeight:
                deg[v1] += edgeToWeight[v2,v1]
    
    D=[]
    A=[]
    for v1 in V:
        d_array = []
        a_array = []
        for v2 in V:
            if v1==v2:
                d_array.append(deg[v1])
                a_array.append(0)
            else:
                d_array.append(0)
                if (v1,v2) in edgeToWeight:
                    a_array.append(edgeToWeight[v1,v2])
                elif (v2,v1) in edgeToWeight:
                    a_array.append(edgeToWeight[v2,v1])
                else:
                    a_array.append(0)
        D.append(d_array)
        A.append(a_array)

    array2d =[]
    for v1 in V:
        array = []
        for v2 in V:
            if v1==v2:
                array.append( deg[v1] )
            elif (v1, v2) in edgeToWeight:
                array.append( -edgeToWeight[v1,v2] )
            elif (v2,v1) in edgeToWeight:
                array.append( -edgeToWeight[v2,v1] )
            else:
                array.append( 0 )
        array2d.append( array )

    PrintMatrix(D)
    print()
    PrintMatrix(A)
    print()
    PrintMatrix(array2d)

    return numpy.array(array2d)


V = [i for i in range(10)]
edge_lists=[
    [1,2,5,8,9],
    [2],
    [],
    [4,5],
    [5],
    [6,7],
    [7],
    [],
    [9],
    []
]
edges = []
for i, l in enumerate(edge_lists):
    for j in l:
        edges.append((i,j))

edgeToWeight={ e:1 for e in edges}

L_nup = GraphLaplacian( V, edgeToWeight )

eig_vals,eig_vecs = numpy.linalg.linalg.eig(L_nup)

eig_vecs = eig_vecs[:,numpy.argsort(eig_vals)]
eig_vals = eig_vals[numpy.argsort(eig_vals)]

for x in list(eig_vals):
    print(f"{x:.3}, ")
for l in list(eig_vecs[:,1:4]):
    for i, x in enumerate(l):
        print(f"{x:.3}",end="")
        if i==len(l)-1:
            print("\\\\")
        else:
            print(" & ",end="")

from matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

ax = fig.add_subplot( projection='3d' )
ax.scatter(eig_vecs[:,1], eig_vecs[:,2],eig_vecs[:,3])

plt.savefig("3dscatter.png")