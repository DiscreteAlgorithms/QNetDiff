def Step0_CalcRelativeAbundance( bacterias, stages, bactToStageToCnts ):
    bactToStageToAbds = dict()
    for b in bacterias:
        bactToStageToAbds[b] = dict()
        for s in stages:
            bactToStageToAbds[b][s] = [ 0 for _ in range( len(bactToStageToCnts[b][s]) ) ]
            
    for s in stages:
        for i in range( len( bactToStageToCnts[bacterias[0]][s] ) ):
            sum = 0.0
            for b in bacterias:
                sum += bactToStageToCnts[b][s][i]
            for b in bacterias:
                bactToStageToAbds[b][s][i] = bactToStageToCnts[b][s][i] / sum
    return bactToStageToAbds