from utility import RelativePathToAbs

def Step1(EDGE_THRESHOLD, stages, bacterias):
    stageToEdgeToCorr = dict()
    for s in stages:
        # read SparCC output
        with open( RelativePathToAbs(f'Correlation/SparCCout_{s}.txt'), "r" )as corrFile:
        #    open(f"for_houkoku_output/correlation.txt","w") as for_texF
            corrFileLines = corrFile.readlines()
            first = True
            stageToEdgeToCorr[s] = dict()
            cnt = 0
            for l in corrFileLines:
                words = l.split()
                if first:
                    # check ここで一行目とbacteriasを照合する
                    first = False
                    # for b in bacterias[:5]:
                    #     for_texF.write(f"&{b}")
                    # for_texF.write(f"\\\\ \\hline\n")
                    continue
                b1 = words[0]
                for i in range( len(bacterias) ):
                    corr = float(words[i+1])
                    b2 = bacterias[i]
                    stageToEdgeToCorr[s][b1, b2] = corr
                # if cnt < 10:
                #     for_texF.write(f"{words[0]}")
                #     for w in words[1:6]:
                #         for_texF.write(f"&{float(w):.3f}")
                #     for_texF.write(f"\\\\ \\hline\n")
                cnt+=1

    adjacency_matrix = dict()
    having_edge = dict()
    for b in bacterias:
        having_edge[b] = False
    for s in stages:
        adjacency_matrix[s] = dict()
        for b1 in bacterias:
            for b2 in bacterias:
                if b1 == b2:
                    continue
                if stageToEdgeToCorr[s][b1, b2] > EDGE_THRESHOLD:
                    adjacency_matrix[s][b1,b2] = stageToEdgeToCorr[s][b1, b2]
                    having_edge[b1] = True
    bacteria_effective = []
    for b in bacterias:
        if having_edge[b]:
            bacteria_effective.append(b)
    return adjacency_matrix, bacteria_effective, stageToEdgeToCorr

def Output_AdjacencyMatrix(EDGE_THRESHOLD, stages, adjacency_matrix, bacterias):
    for s in stages:
        str = ""
        cnt = 0
        for b1 in bacterias:
            for b2 in bacterias:
                if b1 >= b2:
                    continue
                edge_weight = adjacency_matrix[s][b1,b2]
                if adjacency_matrix[s][b1,b2] > EDGE_THRESHOLD:
                    cnt +=1
                    str += f"{b1} {b2} {edge_weight}\n"
        with open(f"tool_output/adjacency_matrix_{s}.txt","w") as adj_file:
            adj_file.write(f"{len(bacterias)} vertices, {cnt} edges\n")
            adj_file.write(str)
