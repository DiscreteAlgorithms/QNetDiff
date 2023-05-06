from utility import RelativePathToAbs
import os.path

def ReadCategory( categoryFileName ):
    sup_category = dict()
    with open( categoryFileName, "r") as file:
        for l in file.readlines():
            words = l.split()
            bacteria = words[0]
            category = words[1]
            sup_category[ bacteria ] = category
    return sup_category

def ReadCntTable( tableFilePath ):
    tableFileName = os.path.splitext( os.path.basename(tableFilePath) )[0]
    stageName = tableFileName
    with open( tableFilePath, "r" ) as tableFile:
        lines = tableFile.readlines()
        
        bactToStageToCnts = dict()
        bacterias = []
        for l in lines[1:]: # skip first line because it has only sample name
            words = l.split()
            bact = words[0]
            if bact not in bacterias:
                bacterias.append( bact )
            if bact not in bactToStageToCnts[ bact ]:
                bactToStageToCnts[ bact ] = dict()
            bactToStageToCnts[ bact ][ stageName ] = dict()
            for w in words[1:]:
                cnt =float( w )
                bactToStageToCnts[ bact ][ stageName ].append( cnt )
    
    return bactToStageToCnts, stageName
            

def ReadCntTable_Whole( tableFilePath, sampleFilePath ):
    tableFile  = open( RelativePathToAbs(tableFilePath), "r" )
    sampleFile = open( RelativePathToAbs(sampleFilePath), "r" )

    sampleFileLines = sampleFile.readlines()
    sampleToStage = dict()
    for l in sampleFileLines:
        sampleData = l.split()
        sampleName = sampleData[0]
        sampleStage = sampleData[1]
        sampleToStage[sampleName] = sampleStage
    
    tableFileLines = tableFile.readlines()
    samples = tableFileLines[0].split()

    bactToStageToCnts = dict()
    bacterias = []
    count_sum = dict()
    for l in tableFileLines[1:len(tableFileLines)]:
        words = l.split()
        bactName = words[0]
        bacterias.append(bactName)
        bactToStageToCnts[bactName] = dict()
        for i in range( 1, len(words) ):
            stage = sampleToStage[samples[i-1]]
            count = float( words[i] )
            if stage not in bactToStageToCnts[bactName].keys():
                bactToStageToCnts[bactName][stage] = []
            bactToStageToCnts[bactName][stage].append(count)
    stages = list( bactToStageToCnts[bacterias[0]].keys() )
    
    stageToSamples = dict()
    for sample in samples:
        stage = sampleToStage[sample]
        if stage not in stageToSamples.keys():
            stageToSamples[stage] = []
        stageToSamples[stage].append(sample)

    return bacterias, stages, bactToStageToCnts, stageToSamples
