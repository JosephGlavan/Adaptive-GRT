# Test simulation to validate conversion from R to Python

# To Do:
# - import stuff (e.g., psychopy, new handler file)
# - write subject model
# - write adaptive block
# - write main exp block
# - data output

import AGRT
import os
from psychopy import data

nReps = 100

#def sm (stimulus):
#    return AGRT.GRTSubjectModel (stimulus, )


# create and open file
if not os.path.isdir('data'):
    os.makedirs('data') #if this fails (e.g. permissions) we will get error
filename = 'data' + os.path.sep + data.getDateStr()
file = filename+'.txt'
logFile=open(file, 'a')

# write header
logFile.write("")

def apf (stimulus, optArgs, result):
    # log stuff
    pass

def gpf (stimulus, optArgs, result):
    # log stuff
    pass

for i in range(nReps):
    AGRT.RunAdaptiveGRTExperiment(trialFunction=AGRT.GRTSubjectModel, nAdaptiveTrials=100, nGRTtrials=1000, 
                                  dim1range=(100, 800), dim2range=(1, 12), dim1steps=100, dim2steps=100, 
                                  lapse=0.0, overallAccuracy=0.75, 
                                  adaptiveFunArgs=((600, 4), (0, 0), (100, .8), (.1, 0), (0, .4), (.9, 0, -.9, 0)), adaptivePostFun=apf, 
                                  blockingFactor=2, 
                                  grtFunArgs=((600, 4), (0, 0), (100, .8), (.1, 0), (0, .4), (.9, 0, -.9, 0)), grtPostFun=gpf)

