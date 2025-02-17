"""
    Adaptive General Recognition Theory Validation Study
    Author: Joseph J Glavan     j.glavan4@gmail.com
    10/28/2021
    
    Two experiments using 2x2 complete identification paradigm with/without adaptive procedure
        Exp1 - Separable dimensions (size and orientation)
        Exp2 - Integral dimensions (saturation and brightness)
    
    
"""

from __future__ import division #so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, monitors, event, logging, gui
from psychopy.constants import * #things like STARTED, FINISHED
from time import sleep
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle, choice
from random import uniform
import math
import os #handy system and path functions
import AGRT

# To Do:
#   - Verify everything that needs to be logged is being logged
#   - Update all software and verify code still works
#   - Is HSV color displaying correctly?
#   - Full test


##########################################################################################################

expVersion = '1.1.0'

##########################################################################################################

# Note: To adjust the size, orientation, saturation, and brightness use the _SIZE, _ORI, STIM_SAT, and STIM_VAL parameters

NUM_BLOCKS = 125 # 1000 test trials
NUM_PRACTICE_BLOCKS = 20 # 160 practice trials
TRIALS_PER_BLOCK = 8 # will be forced to multiples of UNIQUE_STIM_PER_BLOCK if not already
UNIQUE_STIM_PER_BLOCK = 4
NUM_TRIALS = NUM_BLOCKS*TRIALS_PER_BLOCK
BREAK_TRIALS = int(NUM_TRIALS/2)
MAX_TIME = 10   # Maximum duration of a trial (in seconds)
DISPLAY_DIMENSIONS = (1280, 800) #monitors.Monitor('monitor').getSizePix() # Isn't working
ONSET_RANGE = (.25, .75) #(in seconds)
_SIZE = (55, 70) # Features for dimension 1 (sizes in pixels)
_ORI = (50, 70) # Features for dimension 2 (I believe pixels per cycle)
STIM_SIZES = (_SIZE[0], _SIZE[0], _SIZE[1], _SIZE[1]) # Radius of circle stimuli
STIM_ORI = (_ORI[0], _ORI[1], _ORI[0], _ORI[1])
#GAUSS_SD = 10 # Standard deviation of the Gaussian mask for the stimuli
STIM_COL = [-1,-1,-1]
MASK_DURATION = .25
PROBE_DURATION = .25
STIM_LW = 3 # Line width of stimuli
STIM_HUE = 180
STIM_SAT = (.5, .6)
STIM_VAL = (.4, .5)
STIM_COLS = [[STIM_HUE, STIM_SAT[0], STIM_VAL[0]], 
    [STIM_HUE, STIM_SAT[0], STIM_VAL[1]], 
    [STIM_HUE, STIM_SAT[1], STIM_VAL[0]], 
    [STIM_HUE, STIM_SAT[1], STIM_VAL[1]]]
STIM_SIZE = 200
MASK_SIZE = 256
INSTRUCTION_HEIGHT = 40 # size of the font for instructions
ACCURACY_CRITERION = 0.75



#store info about the experiment session
expName='AGRT'    # Name of the Experiment
expInfo={'participant':'TEST', 'exp':['s','i'], 'cond':['a','p']}
dlg=gui.DlgFromDict(dictionary=expInfo,title=expName)
if dlg.OK==False: core.quit() #user pressed cancel
expInfo['date']=data.getDateStr()#add a simple timestamp
expInfo['expName']=expName

#setup files for saving
if not os.path.isdir('data'):
    os.makedirs('data') #if this fails (e.g. permissions) we will get error
filename = 'data' + os.path.sep + '%s_%s' %(expInfo['participant'], expInfo['date'])
file = filename+'.txt'
logFile=open(file, 'a') 

#an ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version=expVersion,
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=False, saveWideText=False,
    dataFileName=filename)

# Setup the Window
win = visual.Window(size=DISPLAY_DIMENSIONS, fullscr=True, screen=0, allowGUI=False, allowStencil=False, color=[0,0,0], colorSpace='rgb', units="pix", monitor="monitor")
# Initialise components for Experiment
globalClock=core.Clock() #to track the time since experiment started
instructionsClock=core.Clock()
trialClock=core.Clock()
responseClock=core.Clock()

def toggle (fc, mode):
    """ Calls .setAutoDraw(mode) for each component in fc. Useful for conglomerate stimuli."""
    for component in fc:
        component.setAutoDraw(mode)

class OrientedCircle (object):
    def __init__(self, win, orientation, rad, col, x, y):
        self._circleStim = visual.Circle(win=win, radius=rad, lineColor=col, fillColor=win.color, pos=[x,y], units='pix', lineColorSpace='rgb', fillColorSpace='rgb', opacity=1.0, lineWidth=STIM_LW)
        #self._lineStim = visual.Line(win=win, lineWidth=STIM_LW, lineColor=col, start=(x,y-rad), end=(x,y+rad), ori=orientation, lineColorSpace='rgb255')
        self._lineStim = visual.Line(win=win, lineWidth=STIM_LW, lineColor=col, start=(x + rad * cos(deg2rad(orientation)), y + rad * sin(deg2rad(orientation))), end=(x + rad * cos(deg2rad(orientation+180)), y + rad * sin(deg2rad(orientation+180))), lineColorSpace='rgb')
        self.pos = [x,y]
        self.radius = rad
        self.color = col
        self.ori = orientation
        self.status = None
    
    def draw(self):
        self._circleStim.draw()
        self._lineStim.draw()
    def setAutoDraw(self, bool):
        self._circleStim.setAutoDraw(bool)
        self._lineStim.setAutoDraw(bool)
    def contains(self, x, y=None):
        if y is None:
            return self._circleStim.contains(x)
        else:
            return self._circleStim.contains(x,y)
    def setPos(self, x, y=None):
        if y is None:
            self._circleStim.setPos(x)
            #self._lineStim.setPos(x)
            self._lineStim.setStart((x[0] + self.radius * cos(deg2rad(self.ori)), x[1] + self.radius * sin(deg2rad(self.ori))))
            self._lineStim.setEnd((x[0] + self.radius * cos(deg2rad(self.ori+180)), x[1] + self.radius * sin(deg2rad(self.ori+180))))
            self.pos = x
        else:
            self._circleStim.setPos(x,y)
            #self._lineStim.setPos(x,y)
            self._lineStim.setStart((x + self.radius * cos(deg2rad(self.ori)), y + self.radius * sin(deg2rad(self.ori))))
            self._lineStim.setEnd((x + self.radius * cos(deg2rad(self.ori+180)), y + self.radius * sin(deg2rad(self.ori+180))))
            self.pos = [x,y]
        

separableStimuli = [OrientedCircle(win, STIM_ORI[0], STIM_SIZES[0], STIM_COL, 0, 0),
    OrientedCircle(win, STIM_ORI[1], STIM_SIZES[1], STIM_COL, 0, 0),
    OrientedCircle(win, STIM_ORI[2], STIM_SIZES[2], STIM_COL, 0, 0),
    OrientedCircle(win, STIM_ORI[3], STIM_SIZES[3], STIM_COL, 0, 0)]
integralStimuli = [visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[0], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
    visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[1], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
    visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[2], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
    visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[3], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv')]

    
# Not sure why this was breaking everything else...
#def findMaskSize (x, scalar):
#   """Takes an iterable x, scales the max value in x by scalar, then finds the smallest square power of two that will cover it."""
#   return (2**math.ceil(math.log(math.sqrt(math.ceil(scalar*max(x))),2)))**2
def newMask ():
#   mysize = findMaskSize(_SIZE, 1.5)
    return visual.GratingStim(win=win, tex=np.random.rand(1024,1024), size=(MASK_SIZE, MASK_SIZE), pos=(0,0), name="mask")
mask1status = "NOT STARTED"; mask2status = "NOT STARTED"


outcomeMessage = visual.TextStim(win=win, ori=0, name='Outcome_message',
    text='You have finished the experiment. Please go see the experimenter at this time.\n\nThank you for your participation!\n\n\n\n\n(Experimenter: Press the spacebar to close)',
    font='Arial', pos=[0, 0], height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, units="pix",
    color='black', colorSpace='rgb', opacity=1, depth=0.0, anchorHoriz='center')

feedbackStimuli = [visual.TextStim(win=win, text="f", units="pix", font='Arial', pos=[0,0], height=200, color="green"),
    visual.TextStim(win=win, text="e", units="pix", font='Arial', pos=[0,0], height=200, color="green"),
    visual.TextStim(win=win, text="j", units="pix", font='Arial', pos=[0,0], height=200, color="green"),
    visual.TextStim(win=win, text="i", units="pix", font='Arial', pos=[0,0], height=200, color="green")]

mouse = event.Mouse(visible=False, win=win)


##### STIMULI TESTING #####
"""
if expInfo['exp'] == 's':
    stimuli = separableStimuli
elif expInfo['exp'] == 'i':
    stimuli = integralStimuli
else:
    raise RuntimeError("Experiment not recognized")
for stim in stimuli:
    stim.draw()
    win.flip()
    while not event.getKeys(["space"]):
        if event.getKeys(["escape"]):
            core.quit()
core.quit()
"""
##### END STIMULI TESTING #####


logFile.write("Version\tDate\tSubject\tExperiment\tCondition\tBlock\tTrial\tDim1\tDim2\tStimulus\tStimNum\tKey\tResponse\tRT\tCorrect\tLambda\tStimEst\n")


def Instructions ():
    initialInstructions.draw()
    win.flip()
    while event.getKeys('space') == []:
        if event.getKeys(["escape"]):
            core.quit()
    for i in range(4):
        exampleInstructions[i].draw()
        stimuli[i].draw()
        win.flip()
        while True:
            theseKeys = event.getKeys()
            if ('f','e','j','i')[i] in theseKeys:
                break
    finalInstructions.draw()
    win.flip()
    while event.getKeys('space') == []:
        if event.getKeys(["escape"]):
            core.quit()

def RunTrial(stimNum, feedback):
    """Runs a single trial with/out feedback and returns the response number"""
    # If we're using a single adaptive block (adapting both dimensions simultaneously) then it needs to return the response
    # If we're using two adaptive blocks (adapting to each dimension separately) then it only needs to return correct/incorrect
    global stimuli
    
    stimuli[stimNum].status = "NOT STARTED"
    stimuli[stimNum].setAutoDraw(False)
    #mask1status = "NOT STARTED"
    mask2status = "NOT STARTED"
    response = None; rt = -1
    RSOI = uniform(ONSET_RANGE[0], ONSET_RANGE[1])
    t = 0; trialClock.reset(); continueRoutine = True
    while continueRoutine:
        t = trialClock.getTime()
        if t > MAX_TIME:
            continueRoutine = False
        if t >= RSOI and stimuli[stimNum].status == "NOT STARTED":
            stimuli[stimNum].setAutoDraw(True)
            stimuli[stimNum].status = "STARTED"
            responseClock.reset()
            event.clearEvents()
        if t >= RSOI + PROBE_DURATION and stimuli[stimNum].status == "STARTED":
            stimuli[stimNum].setAutoDraw(False)
            stimuli[stimNum].status = "FINISHED"
        if t >= RSOI + PROBE_DURATION and mask2status=="NOT STARTED":
            mask = newMask()
            mask.setAutoDraw(True)
            mask2status = "STARTED"
        if t >= RSOI + PROBE_DURATION + MASK_DURATION and mask2status=="STARTED":
            mask.setAutoDraw(False)
            mask2status = "FINISHED"
        if t >= RSOI and stimuli[stimNum].status != "NOT STARTED":
            keys = event.getKeys(keyList=['f','e','j','i'], timeStamped=responseClock)
            for k in keys:
                if rt < k[1]:
                    response = k[0]
                    rt = k[1]
        if response is not None and mask2status == "FINISHED":
            if feedback == True:
                continueRoutine = False
            elif t >= RSOI + PROBE_DURATION + MASK_DURATION + 1:
                continueRoutine = False
                
        if event.getKeys(["escape"]):
            core.quit()
        win.flip()
            
    # Score the Trial
    respnum = -1
    if response == 'f':
        respnum = 0
    elif response == 'e':
        respnum = 1
    elif response == 'j':
        respnum = 2
    elif response == 'i':
        respnum = 3
    if respnum == stimNum:
        mycorrect = 1
    else:
        mycorrect = 0
        
    # Provide Feedback
    if feedback == True:
        if mycorrect == 1:
            feedbackStimuli[stimNum].setColor("green")
        else:
            feedbackStimuli[stimNum].setColor("red")
        feedbackStimuli[stimNum].draw()
        trialClock.reset()
        win.flip()
        while trialClock.getTime() < 1:
            if event.getKeys(["escape"]):
                core.quit()
        win.flip()
    
    if expInfo['exp']=='s':
        dim1 = stimuli[stimNum].radius
        dim2 = stimuli[stimNum].ori
    elif expInfo['exp']=='i':
        dim1 = stimuli[stimNum].fillColor[1]
        dim2 = stimuli[stimNum].fillColor[2]
    else:
        raise RuntimeError("Experiment not recognized")
    
    # Log the Trial    
    
    logstr = "\t".join(map(str, [expVersion, expInfo['date'], expInfo['participant'], expInfo['exp'], expInfo['cond'], blockName, trialCounter, dim1, dim2, ((0,0), (0,1), (1,0), (1,1))[stimNum], stimNum, response, respnum, rt, mycorrect, "NA", "NA"])) + "\n"
    logFile.write(logstr)
    #core.wait(1.0)
    return respnum


def RunAgrtTrial (stimulus, logfile=None, handler=None, trial=None, info=None, additionalArgs=None):
    """Modifying RunTrial to work with adaptive code blocks.
    To Do: Integrate all code more cleanly in v2.0.0
    """
    
    # Unpack variables
    stimX,stimY = stimulus
    feedback = additionalArgs[0]
    
    # Check for break
    if trial % 100 == 1 and additionalArgs[-1] == 'main': # Offer a break every 100 trials during main block
        visual.TextStim(win=win, ori=0, name='breakInstructions',text="You may now take a short break.\n\nPress the spacebar when you are ready to continue.",
        font='Arial', pos=[0, 0], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center', color='black', colorSpace='rgb', opacity=1, depth=0.0).draw()
        win.flip()
        while event.getKeys('space') == []:
            if event.getKeys(["escape"]):
                core.quit()
        win.flip()
    
    # Construct Stimulus
    if info['exp'] == 's':
        myStim = OrientedCircle(win, stimY, stimX, STIM_COL, 0, 0)
    elif info['exp'] == 'i':
        myStim = visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=[STIM_HUE, stimX, stimY], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv')
    else:
        raise RuntimeError("Experiment not recognized")
    
    
    myStim.status = "NOT STARTED"
    myStim.setAutoDraw(False)
    #mask1status = "NOT STARTED"
    mask2status = "NOT STARTED"
    response = None; rt = -1
    RSOI = uniform(ONSET_RANGE[0], ONSET_RANGE[1])
    t = 0; trialClock.reset(); continueRoutine = True
    while continueRoutine:
        t = trialClock.getTime()
        if t > MAX_TIME:
            continueRoutine = False
        if t >= RSOI and myStim.status == "NOT STARTED":
            myStim.setAutoDraw(True)
            myStim.status = "STARTED"
            responseClock.reset()
            event.clearEvents()
        if t >= RSOI + PROBE_DURATION and myStim.status == "STARTED":
            myStim.setAutoDraw(False)
            myStim.status = "FINISHED"
        if t >= RSOI + PROBE_DURATION and mask2status=="NOT STARTED":
            mask = newMask()
            mask.setAutoDraw(True)
            mask2status = "STARTED"
        if t >= RSOI + PROBE_DURATION + MASK_DURATION and mask2status=="STARTED":
            mask.setAutoDraw(False)
            mask2status = "FINISHED"
        if t >= RSOI and myStim.status != "NOT STARTED":
            keys = event.getKeys(keyList=['f','e','j','i'], timeStamped=responseClock)
            for k in keys:
                if rt < k[1]:
                    response = k[0]
                    rt = k[1]
        if response is not None and mask2status == "FINISHED":
            if feedback == True:
                continueRoutine = False
            elif t >= RSOI + PROBE_DURATION + MASK_DURATION + 1:
                continueRoutine = False
                
        if event.getKeys(["escape"]):
            core.quit()
        win.flip()
            
    # Score the Trial
    respnum = -1
    responseTuple = None
    if response == 'f':
        respnum = 0
        responseTuple = (0,0)
    elif response == 'e':
        respnum = 1
        responseTuple = (0,1)
    elif response == 'j':
        respnum = 2
        responseTuple = (1,0)
    elif response == 'i':
        respnum = 3
        responseTuple = (1,1)
    
    if isinstance(handler, AGRT.GRTHandler):
        if responseTuple == eval(handler.thisTrial):
            mycorrect = 1
        else:
            mycorrect = 0    
        if eval(handler.thisTrial) == (0,0):
            stimNum = 0
        elif eval(handler.thisTrial) == (0,1):
            stimNum = 1
        elif eval(handler.thisTrial) == (1,0):
            stimNum = 2
        elif eval(handler.thisTrial) == (1,1):
            stimNum = 3
        else:
            raise RuntimeError("stimNum could not be determined.")
    else:
        mycorrect = -1
        stimNum = -1
        
    # Provide Feedback
    if feedback == True:
        if mycorrect == 1:
            feedbackStimuli[stimNum].setColor("green")
        else:
            feedbackStimuli[stimNum].setColor("red")
        feedbackStimuli[stimNum].draw()
        trialClock.reset()
        win.flip()
        while trialClock.getTime() < 1:
            if event.getKeys(["escape"]):
                core.quit()
        win.flip()

    
    # Log the Trial    
        
    logstr = "\t".join(map(str, [expVersion, expInfo['date'], expInfo['participant'], expInfo['exp'], expInfo['cond'], additionalArgs[-1], trial, stimX, stimY, 
                                 handler.thisTrial if isinstance(handler, AGRT.GRTHandler) else "NA", 
                                 stimNum, response, respnum, rt, mycorrect, 
                                 handler.estimateLambda() if isinstance(handler, AGRT.AGRTHandler) else "NA",
                                 handler.estimateGRTintensities(ACCURACY_CRITERION) if isinstance(handler, AGRT.AGRTHandler) else "NA"])) + "\n"
    logFile.write(logstr)
    #core.wait(1.0)
    
    return responseTuple


# Initialize condition parameters
if expInfo['exp'] == 's':
    stimuli = separableStimuli
    instrArgs = ["circle with a line through it", "circle", "smaller", "larger", "less steep", "more steep"]
elif expInfo['exp'] == 'i':
    stimuli = integralStimuli
    instrArgs = ["blue rectangle", "rectangle", "less saturated", "more saturated", "less bright", "more bright"]
else:
    raise RuntimeError("Experiment not recognized")
if expInfo['cond'] == 'a':
    # Initialize the PsiHandler
    instrArgs.append("without feedback")
elif expInfo['cond'] == 'p':
    instrArgs.append("with feedback")
else:
    raise RuntimeError("Condition not recognized")
    
    
initialInstructions = visual.TextStim(win=win, ori=0, name='TaskInstructions',
    text="You will be presented with a {obj} followed by a mask.\n\n\n\nPress 'f' if the {obj1} is {dim1a} and {dim2a}.\n\nPress 'e' if the {obj1} is {dim1a} and {dim2b}.\n\nPress 'j' if the {obj1} is {dim1b} and {dim2a}.\n\nPress 'i' if the {obj1} is {dim1b} and {dim2b}.\n\n\n\nWhen you are ready, press the spacebar now to see examples.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5]),
    font='Arial', pos=[0, 0], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0)#, anchorHoriz='left')
exampleInstructions = [visual.TextStim(win=win, ori=0, name='smallLowInstructions',
    text="Press 'f' if the {obj1} is {dim1a} and {dim2a}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5]),
    font='Arial', pos=[0, DISPLAY_DIMENSIONS[1]*.3], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0),
    visual.TextStim(win=win, ori=0, name='smallHighInstructions',
    text="Press 'e' if the {obj1} is {dim1a} and {dim2b}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5]),
    font='Arial', pos=[0, DISPLAY_DIMENSIONS[1]*.3], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0),
    visual.TextStim(win=win, ori=0, name='largeLowInstructions',
    text="Press 'j' if the {obj1} is {dim1b} and {dim2a}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5]),
    font='Arial', pos=[0, DISPLAY_DIMENSIONS[1]*.3], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0),
    visual.TextStim(win=win, ori=0, name='largeHighInstructions',
    text="Press 'i' if the {obj1} is {dim1b} and {dim2b}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5]),
    font='Arial', pos=[0, DISPLAY_DIMENSIONS[1]*.3], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0)]
finalInstructions = visual.TextStim(win=win, ori=0, name='InstructionsReminder',
    text="You will now complete a relatively short practice block {optfb}. Press the spacebar when you are ready to begin the experiment.\n\nREMEMBER:\n\nPress 'f' if the {obj1} is {dim1a} and {dim2a}.\n\nPress 'e' if the {obj1} is {dim1a} and {dim2b}.\n\nPress 'j' if the {obj1} is {dim1b} and {dim2a}.\n\nPress 'i' if the {obj1} is {dim1b} and {dim2b}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5], optfb=instrArgs[6]),
    font='Arial', pos=[0, 0], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
    color='black', colorSpace='rgb', opacity=1, depth=0.0)#, anchorHoriz='left')

#### PRACTICE BLOCK ####
blockName = "PRACT"

trialOrder = list(range(UNIQUE_STIM_PER_BLOCK))*(TRIALS_PER_BLOCK//UNIQUE_STIM_PER_BLOCK)

trialCounter = 0
if expInfo['cond'] == 'a':
    # Adaptive Condition
    # Set up extreme stimuli, run 3 blocks (12 trials) with them
    # Run the remaining trials passing responses to the PsiHandler and resetting stimuli between trials
    # Save stimuli
    
    
    # Override 'p' condition stimuli parameters to become the extremes for each dimension
    _SIZE = (20, 100) # Features for dimension 1 (sizes in pixels)
    _ORI = (50, 85) # Features for dimension 2 (I believe pixels per cycle)
    STIM_SIZES = (_SIZE[0], _SIZE[0], _SIZE[1], _SIZE[1]) # Radius of circle stimuli
    STIM_ORI = (_ORI[0], _ORI[1], _ORI[0], _ORI[1])
    STIM_SAT = (.3, .7)
    STIM_VAL = (.3, .7)
    STIM_COLS = [[STIM_HUE, STIM_SAT[0], STIM_VAL[0]], 
        [STIM_HUE, STIM_SAT[0], STIM_VAL[1]], 
        [STIM_HUE, STIM_SAT[1], STIM_VAL[0]], 
        [STIM_HUE, STIM_SAT[1], STIM_VAL[1]]]
    separableStimuli = [OrientedCircle(win, STIM_ORI[0], STIM_SIZES[0], STIM_COL, 0, 0),
        OrientedCircle(win, STIM_ORI[1], STIM_SIZES[1], STIM_COL, 0, 0),
        OrientedCircle(win, STIM_ORI[2], STIM_SIZES[2], STIM_COL, 0, 0),
        OrientedCircle(win, STIM_ORI[3], STIM_SIZES[3], STIM_COL, 0, 0)]
    integralStimuli = [visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[0], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
        visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[1], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
        visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[2], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv'),
        visual.Rect(win, units='pix', lineWidth=0, lineColor=None, fillColor=STIM_COLS[3], pos=(0, 0), size=STIM_SIZE, fillColorSpace='hsv')]
    
    if expInfo['exp'] == 's':
        stimuli = separableStimuli
        dim1Extremes = _SIZE
        dim2Extremes = _ORI
    elif expInfo['exp'] == 'i':
        stimuli = integralStimuli
        dim1Extremes = STIM_SAT
        dim2Extremes = STIM_VAL
    else:
        raise RuntimeError("Experiment not recognized")
    
    Instructions()
    for block in range(2):
        shuffle(trialOrder)
        for trial in trialOrder:
            trialCounter += 1
            RunTrial(trial, True)

    
    AGRT.RunAdaptiveGRTExperiment(trialFunction = RunAgrtTrial, nAdaptiveTrials = 8, nGRTtrials = 8, #144, nGRTtrials = NUM_TRIALS, 
                                  dim1range = dim1Extremes, dim2range = dim2Extremes, 
                                  dim1steps=100, dim2steps=100, lapse=0.0, 
                                  overallAccuracy=ACCURACY_CRITERION, blockingFactor=2, 
                                  adaptiveFunArgs=[False], 
                                  grtFunArgs=[False],
                                  info=expInfo, logfile=logFile, savePosterior=True)
    
elif expInfo['cond'] == 'p':
    # Pilot Condition
    Instructions()
    
    for block in range(NUM_PRACTICE_BLOCKS):
        shuffle(trialOrder)
        for trial in trialOrder:
            trialCounter += 1
            RunTrial(trial, True)
    
    #### MAIN BLOCK ####
    blockName = "MAIN"
    
    visual.TextStim(win=win, ori=0, name='InstructionsReminder',
        text="You will now complete the main experiment without feedback. You will periodically be offered a break. Press the spacebar when you are ready to begin.\n\nREMEMBER:\n\nPress 'f' if the {obj1} is {dim1a} and {dim2a}.\n\nPress 'e' if the {obj1} is {dim1a} and {dim2b}.\n\nPress 'j' if the {obj1} is {dim1b} and {dim2a}.\n\nPress 'i' if the {obj1} is {dim1b} and {dim2b}.".format(obj=instrArgs[0], obj1=instrArgs[1], dim1a=instrArgs[2], dim1b=instrArgs[3], dim2a=instrArgs[4], dim2b=instrArgs[5], optfb=instrArgs[6]),
        font='Arial', pos=[0, 0], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center',
        color='black', colorSpace='rgb', opacity=1, depth=0.0).draw()
    win.flip()
    while event.getKeys('space') == []:
        if event.getKeys(["escape"]):
            core.quit()
    win.flip()
    
    
    trialCounter = 0
    for block in range(NUM_BLOCKS):
    
        if block % 25 == 0 and block != 0 and block != NUM_TRIALS: # Offer a break every 100 trials
            visual.TextStim(win=win, ori=0, name='breakInstructions',text="You may now take a short break.\n\nPress the spacebar when you are ready to continue.",
            font='Arial', pos=[0, 0], units="pix", height=INSTRUCTION_HEIGHT,wrapWidth=DISPLAY_DIMENSIONS[0]-100, anchorVert='center', color='black', colorSpace='rgb', opacity=1, depth=0.0).draw()
            win.flip()
            while event.getKeys('space') == []:
                if event.getKeys(["escape"]):
                    core.quit()
            win.flip()
    
        shuffle(trialOrder)
        for trial in trialOrder:
            trialCounter += 1
            RunTrial(trial, False)

else:
    raise RuntimeError("Condition not recognized")

outcomeMessage.draw()
win.flip()
while event.getKeys('space') == []:
    if event.getKeys(["escape"]):
        core.quit()
win.flip()

logFile.close()
win.close()
core.quit()




