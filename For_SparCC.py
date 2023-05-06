import subprocess

def Run_SparCC( corrFileName_withOutStageAndExtension, stages, bacterias, bactToStageToCnts, stageToSamples):
    for s in stages:
        # make input file for SparCC
        fileName = f'SparCC_Input/SparCCInput_{s}.tsv'
        with open( fileName, 'w' ) as file:
            for sampleName in stageToSamples[s]:
                file.write(f"{sampleName}\t" )
            file.write("\n")
            for b in bacterias:
                file.write(b)
                for cnt in bactToStageToCnts[b][s]:
                    # file.write('\t{:.12f}'.format(abd).rstrip('0') )
                    file.write(f"\t{cnt}")
                file.write('\n')
        # run SparCC
        subprocess.run( ['python', 'SparCC3-master/Sparcc.py', fileName, '-c' , f'Correlation/{corrFileName_withOutStageAndExtension}{s}.cor', '-v', f'SparCC_Output/{s}.cov'] )
        print('wrote ' + f'SparCC_Output/SparCCout_{s}.txt')
        print("=================")
  
