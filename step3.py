from scipy import stats

# return top10 by Mann-Whitney U test
def Step3( MWU_P_TH, MWU_RANK_TH, baseStage, focusStage, bactToStageToAbds):
    p_bactVec = []
    bactToP_val = dict()
    minInOver_p_bact = 1,""
    for b in bactToStageToAbds.keys():
        x = bactToStageToAbds[b][focusStage]
        y = bactToStageToAbds[b][baseStage]
        # hypothesis testing
        p_val = stats.mannwhitneyu(x,y,alternative='greater').pvalue
        bactToP_val[b] = p_val
        if p_val > MWU_P_TH:
            minInOver_p_bact = min( minInOver_p_bact, (p_val,b) )
            continue
        p_bactVec.append( (p_val, b) )

    bacteria_focused = [ sorted(p_bactVec)[i][1] for i in range( min( MWU_RANK_TH, len(p_bactVec) ) )]
    
    return bacteria_focused, bactToP_val        

