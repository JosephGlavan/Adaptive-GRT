# Joseph J Glavan 2023-09-15
# Test simulation to validate conversion from R to Python
# Simulates the same subject 100 times through 100 adaptive trials and 1000 GRT trials

import AGRT
import os
from psychopy import data, core


nReps = 100


# Create and open file
if not os.path.isdir('data'):
    os.makedirs('data') #if this fails (e.g. permissions) we will get error
filename = 'data' + os.path.sep + data.getDateStr()
file = filename+'.txt'
logFile=open(file, 'a')

# Write header
logFile.write("Version\tDate\tSubject\tCondition\tParams\tBlock\tTrial\tStimulus\tResponse\tRT\tLambda\n")

# Simulate the study
for i in range(nReps):
    expInfo = [1.0, # Version
               data.getDateStr(), # Date
               i, # Participant
               'a'] # Cond
    
    AGRT.RunAdaptiveGRTExperiment(trialFunction=AGRT.GRTSubjectModel, nAdaptiveTrials=100, nGRTtrials=1000, 
                                  dim1range=(100, 800), dim2range=(1, 12), dim1steps=100, dim2steps=100, 
                                  lapse=0.0, overallAccuracy=0.75, blockingFactor=2, 
                                  adaptiveFunArgs=((600, 4), (0, 0), (100, .8), (.1, 0), (0, .4), (.9, 0, -.9, 0)), 
                                  grtFunArgs=((600, 4), (0, 0), (100, .8), (.1, 0), (0, .4), (.9, 0, -.9, 0)),
                                  info=expInfo, logfile=logFile)

logFile.close()
core.quit()