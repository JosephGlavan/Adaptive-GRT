# =============================================================================
#   Copyright 2022 Joseph J Glavan
# 
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# 
#   Please see <http://www.gnu.org/licenses/> for a copy of the GNU General Public License.
# =============================================================================

from __future__ import absolute_import, division, print_function
from builtins import range
from builtins import object
#import math
import warnings
import random
import os
#import sys
#import time
import numpy as np
#import scipy
from scipy import stats, special
from scipy.stats.mvn import mvnun
from psychopy import logging
from psychopy.tools.fileerrortools import handleFileCollision
from psychopy.data import TrialHandler
from psychopy.data.staircase import StairHandler
__all__ = ['AGRTHandler','RunAdaptiveBlock','GRTHandler','RunGRTBlock','RunAdaptiveGRTExperiment','GRTSubjectModel']


# =============================================================================
# Internal Psi Object for AGRT Handler
# =============================================================================

class agrtPsiObject(object):
    """
    Special class to handle internal array and functions of Psi adaptive psychophysical method (Kontsevich & Tyler, 1999).
    Modified to use Adaptive GRT measurement model.
    """
    
    def __init__(self, x, alpha, beta, xPrecision, aPrecision, bPrecision, delta=0, stepType='lin', prior=None):
        self._probLambdaGivenXResponse = self._probResponseGivenX = self._probLambdaGivenXResponse = self._entropyXResponse = self._expectedEntropyX = self.nextIntensityIndex = self.nextIntensity = None

        #Save dimensions
        if stepType == 'lin':
            self.x = np.linspace(x[0], x[1], xPrecision, True)
        elif stepType == 'log':
            self.x = np.logspace(np.log10(x[0]), np.log10(x[1]), xPrecision, True)
        else:
            raise RuntimeError('Invalid step type. Unable to initialize PsiObject.')
        self.alpha = np.linspace(alpha[0], alpha[1], aPrecision)
        self.beta = np.linspace(beta[0], beta[1], bPrecision)
        #self.alpha = np.linspace(alpha[0], alpha[1], int(round((alpha[1]-alpha[0])/aPrecision)+1), True)
        #self.beta = np.linspace(beta[0], beta[1], int(round((beta[1]-beta[0])/bPrecision)+1), True)
        self.r = np.array(list(range(2)))
        self.delta = delta
        
        # Change x,a,b,r arrays to matrix computation compatible orthogonal 4D arrays
        # ALWAYS use the order for P(r|lambda,x); i.e. [r,a,b,x]
        self._r = self.r.reshape((self.r.size,1,1,1))
        self._alpha = self.alpha.reshape((1,self.alpha.size,1,1))
        self._beta = self.beta.reshape((1,1,self.beta.size,1))
        self._x = self.x.reshape((1,1,1,self.x.size))
        
        #Create P(lambda)
        if prior is None or prior.shape != (1, len(self.alpha),len(self.beta), 1):
            if prior is not None:
                warnings.warn("Prior has incompatible dimensions. Using uniform (1/N) probabilities.")
            self._probLambda = np.ndarray(shape=(1,len(self.alpha),len(self.beta),1))
            self._probLambda.fill(1/(len(self.alpha)*len(self.beta)))
        else:
            if prior.shape == (1, len(self.alpha), len(self.beta), 1):
                self._probLambda = prior
            else:
                self._probLambda = prior.reshape(1, len(self.alpha), len(self.beta), 1)
                
        #Create P(r | lambda, x)
        
        # self._probResponseGivenLambdaX = np.zeros(shape=(len(self.r), len(self.alpha), len(self.beta), len(self.x)))
        #                            
        # for i in range(len(self.alpha)):
        #     for j in range(len(self.beta)):
        #         for k in range(len(self.x)):
        #             # alpha == the decision threshold
        #             # beta == the sd of the perceptual distribution
        #             # x == the mean of the perceptual distribution
        #             self._probResponseGivenLambdaX[0,i,j,k] = stats.norm.cdf(self.alpha[i], loc=self.x[k], scale=self.beta[j])
        #             self._probResponseGivenLambdaX[1,i,j,k] = 1 - stats.norm.cdf(self.alpha[i], loc=self.x[k], scale=self.beta[j])

        # Without lapse
        #self._probResponseGivenLambdaX = np.array([0,1]).reshape(2,1,1,1) + np.array([1,-1]).reshape(2,1,1,1) * stats.norm.cdf(self._alpha, loc=self._x, scale=self._beta)
        # With lapse
        self._probResponseGivenLambdaX = np.array([0,1]).reshape(2,1,1,1) + np.array([1,-1]).reshape(2,1,1,1) * ((self.delta/2) + (1 - self.delta) * stats.norm.cdf(self._alpha, loc=self._x, scale=self._beta))

#        if TwoAFC:
#            self._probResponseGivenLambdaX = (1-self._r) + (2*self._r-1) * ((.5 + .5 * stats.norm.cdf(self._x, self._alpha, self._beta)) * (1 - self.delta) + self.delta / 2)
#        else: # Yes/No
#            self._probResponseGivenLambdaX = (1-self._r) + (2*self._r-1) * (stats.norm.cdf(self._x, self._alpha, self._beta)*(1-self.delta)+self.delta/2)
        
    def update(self, response=None):
        if response is not None:    #response should only be None when Psi is first initialized
            self._probLambda = self._probLambdaGivenXResponse[response,:,:,self.nextIntensityIndex].reshape((1,len(self.alpha),len(self.beta),1))
            
        #Create P(r | x)
        self._probResponseGivenX = np.sum(self._probResponseGivenLambdaX * self._probLambda, axis=(1,2)).reshape((len(self.r),1,1,len(self.x)))
        
        #Create P(lambda | x, r)
        self._probLambdaGivenXResponse = self._probLambda*self._probResponseGivenLambdaX/self._probResponseGivenX
        
        #Create H(x, r)
        self._entropyXResponse = -1* np.sum(self._probLambdaGivenXResponse * np.log10(self._probLambdaGivenXResponse, out = np.zeros_like(self._probLambdaGivenXResponse), where = self._probLambdaGivenXResponse > 0), axis=(1,2)).reshape((len(self.r),1,1,len(self.x)))
        
        #Create E[H(x)]
        self._expectedEntropyX = np.sum(self._entropyXResponse * self._probResponseGivenX, axis=0).reshape((1,1,1,len(self.x)))
        
        #Generate next intensity
        self.nextIntensityIndex = np.argmin(self._expectedEntropyX, axis=3)[0][0][0]
        self.nextIntensity = self.x[self.nextIntensityIndex]
        
    def estimateLambda(self):
        return (np.sum(np.sum(self._alpha.reshape((len(self.alpha),1))*self._probLambda.squeeze(), axis=1)), np.sum(np.sum(self._beta.reshape((1,len(self.beta)))*self._probLambda.squeeze(), axis=1)))
        
    def estimateThreshold(self, thresh, lam):
        if lam is None:
            lamb = self.estimateLambda()
        else:
            lamb = lam

        # With lapse
        return (lamb[0] - lamb[1] * np.sqrt(2) * special.erfinv((2 * (np.sqrt(thresh)) - self.delta) / (1 - self.delta) - 1),
                lamb[0] - lamb[1] * np.sqrt(2) * special.erfinv((2 * (1 - np.sqrt(thresh)) - self.delta) / (1 - self.delta) - 1))
       
#         Without lapse
#         return (lamb[0] - lamb[1] * np.sqrt(2) * special.erfinv(2 * np.sqrt(thresh) - 1),
#                 lamb[0] - lamb[1] * np.sqrt(2) * special.erfinv(2 * (1-np.sqrt(thresh)) - 1))

#        if self._TwoAFC:
#            return stats.norm.ppf((2*thresh-1)/(1-self.delta), lamb[0], lamb[1])
#        else:
#            return stats.norm.ppf((thresh-self.delta/2)/(1-self.delta), lamb[0], lamb[1])
        
    def savePosterior(self): #, file):
        #np.save(file, self._probLambda) # commented out because now saving the actual file in the AGRTHandler - only need the p(lambda) from each psi object
        return self._probLambda


# =============================================================================
# AGRT Handler
# =============================================================================

# Copy/paste the original PsiHandler code from PsychoPy, may have to find an older version
# Edit it to initialize like the R Adapt function etc.

class AGRTHandler (StairHandler):
    """
    Handler object to manage the two internal Psi objects
    """
    
    def __init__(self,
                 nTrials, dim1range, dim2range,
                 #delta1, delta2=None,
                 #dim1stepType='lin', dim2stepType='lin',
                 dim1steps=100, dim2steps=100,
                 lapse=0,
                 prior=None,
                 fromFile=False,
                 extraInfo=None,
                 name=''):
        # UPDATE DOC STRING
        """Initializes the handler and creates an internal Psi Object for
        grid approximation.
        :Parameters:
            nTrials (int)
                The number of trials to run.
            intensRange (list)
                Two element list containing the (inclusive) endpoints of
                the stimuli intensity range.
            alphaRange  (list)
                Two element list containing the (inclusive) endpoints of
                the alpha (location parameter) range.
            betaRange   (list)
                Two element list containing the (inclusive) endpoints of
                the beta (slope parameter) range.
            intensPrecision (float or int)
                If stepType == 'lin', this specifies the step size of the
                stimuli intensity range. If stepType == 'log', this specifies
                the number of steps in the stimuli intensity range.
            alphaPrecision  (float)
                The step size of the alpha (location parameter) range.
            betaPrecision   (float)
                The step size of the beta (slope parameter) range.
            delta   (float)
                The guess rate.
            stepType    (str)
                The type of steps to be used when constructing the stimuli
                intensity range. If 'lin' then evenly spaced steps are used.
                If 'log' then logarithmically spaced steps are used.
                Defaults to 'lin'.
            expectedMin  (float)
                The expected lower asymptote of the psychometric function
                (PMF).
                For a Yes/No task, the PMF usually extends across the
                interval [0, 1]; here, `expectedMin` should be set to `0`.
                For a 2-AFC task, the PMF spreads out across [0.5, 1.0].
                Therefore, `expectedMin` should be set to `0.5` in this
                case, and the 2-AFC psychometric function described above
                going to be is used.
                Currently, only Yes/No and 2-AFC designs are supported.
                Defaults to 0.5, or a 2-AFC task.
            prior   (numpy ndarray or str)
                Optional prior distribution with which to initialize the
                Psi Object. This can either be a numpy ndarray object or
                the path to a numpy binary file (.npy) containing the ndarray.
            fromFile    (str)
                Flag specifying whether prior is a file pathname or not.
            extraInfo   (dict)
                Optional dictionary object used in PsychoPy's built-in
                logging system.
            name    (str)
                Optional name for the PsiHandler used in PsychoPy's built-in
                logging system.
        :Raises:
            NotImplementedError
                If the supplied `minVal` parameter implies an experimental
                design other than Yes/No or 2-AFC.
        """

        StairHandler.__init__(
            self, startVal=None, nTrials=nTrials, extraInfo=extraInfo,
            #stepType=stepType, minVal=intensRange[0], maxVal=intensRange[1], 
            name=name
        )

        self.finished = None

        # Initialize priors
        if prior is not None and fromFile:
            try:
                prior = np.load(prior)
            except IOError:
                logging.warning("The specified pickle file could not be "
                                "read. Using a uniform prior instead.")
                prior = None
        if prior is None:
            prior = [None,None]
     
        # Calculate beta ranges
        # max beta is the standard deviation that yields ( 99% response probability - half the marginal lapse rate i.e. 0.01 less than the saturation) for a stimulus corresponding to the min of dimRange and decision bound at mean(dimRange)
        # min beta is max beta divided by the number of steps
        
        # With lapse
        dim1betaRange = [0, (np.average(dim1range) - dim1range[0]) / (np.sqrt(2) * special.erfinv((2 * (.99-lapse/2) - lapse) / (1 - lapse) - 1))]
        # Without lapse
        # dim1betaRange = [0, (np.average(dim1range) - dim1range[0]) / (np.sqrt(2) * special.erfinv(2 * .99 - 1))]
        dim1betaRange[0] = dim1betaRange[1]/dim1steps
        # With lapse
        dim2betaRange = [0, (np.average(dim2range) - dim2range[0]) / (np.sqrt(2) * special.erfinv((2 * (.99-lapse/2) - lapse) / (1 - lapse) - 1))]
        # Without lapse
        #dim2betaRange = [0, (np.average(dim2range) - dim2range[0]) / (np.sqrt(2) * special.erfinv(2 * .99 - 1))]
        dim2betaRange[0] = dim2betaRange[1]/dim2steps


        # Create Psi objects
        self._psi1 = agrtPsiObject(
            dim1range, dim1range, dim1betaRange, 
            dim1steps, dim1steps, dim1steps, 
            delta=lapse, stepType='lin', prior=prior[0])
        self._psi1.update(None)
        self._psi2 = agrtPsiObject(
            dim2range, dim2range, dim2betaRange, 
            dim2steps, dim2steps, dim2steps, 
            delta=lapse, stepType='lin', prior=prior[1])
        self._psi2.update(None)        

    def addResponse(self, result, intensity=None):
        """Record a response and update the internal Psi object.
        The result and intensity should be tuples of length equal
        to the number of dimensions.
        
        Result is a tuple of integers reflecting the response level (n,m).
        Intensity is a tuple of floats reflecting the stimumulus intensity on each dimension.
        
        Supplying an `intensity` value here
        indicates that you did not use the
        recommended intensity in your last trial and the staircase will
        replace its recorded value with the one you supplied here.
        """

        # add the current data to experiment if possible
        self.data.append(result)
        if self.getExp() is not None:
            # update the experiment handler too
            self.getExp().addData(self.name + ".response", result)
            
        # if needed replace the existing intensity with this custom one # Not yet supported
        if intensity is not None:
            self.intensities.pop()
            self.intensities.append(intensity)
            # if possible, update the psi objects (might not be if the intensity is not a value in the internal array)
            ## To Do
            raise RuntimeWarning("Specifying an alternative intensity is not currently supported.")

        self._psi1.update(result[0])
        self._psi2.update(result[1])

    def __next__(self):
        """Advances to next trial and returns it.
        """
        self._checkFinished()
        if self.finished == False:
            # update pointer for next trial
            self.thisTrialN += 1
            self.intensities.append((self._psi1.nextIntensity,self._psi2.nextIntensity))
            return (self._psi1.nextIntensity, self._psi2.nextIntensity)
        else:
            self._terminate()

    next = __next__  # allows user to call without a loop `val = trials.next()`

    def _checkFinished(self):
        """checks if we are finished.
        Updates attribute: `finished`
        """
        if self.nTrials is not None and len(self.intensities) >= self.nTrials:
            self.finished = True
        else:
            self.finished = False

    def estimateLambda(self):
        """Returns a tuple of lambdas where each is a tuple of (location, slope)
        """
        return (self._psi1.estimateLambda(),self._psi2.estimateLambda())

    def estimateGRTintensities(self, overallAccuracy, lambdas=None):
        """Returns a tuple of tuples ((L1,H1),(L2,H2)) providing the low (L) 
        and high (H) intensities for each dimension (1,2) based on the overall 
        accuracy specified.

        The optional argument 'lambdas' allows thresholds to be estimated
        without having to recompute the maximum likelihood lambda. It should be 
        provided in the same form as returned by the method 'estimateLambda': 
        ((location1,slope1),(location2,slope2))
        """
        
        if lambdas is None:
            lambdas = (None, None)
        if lambdas != (None, None):
            try:
                if len(lambdas) != 2:
                    msg = f"Invalid user-specified lambdas tuple: {lambdas}. New estimates for lambda will be computed."
                    warnings.warn(msg, SyntaxWarning)
                    lambdas = (None, None)
                else:
                    for i in range(len(lambdas)):
                        if len(lambdas[i]) != 2:
                            msg = f"Invalid user-specified lambda pair: {lambdas}. A new estimate of lambda will be computed."
                            warnings.warn(msg, SyntaxWarning)
                            lambdas[i] = None
            except TypeError:
                msg = f"Invalid user-specified lambdas tuple: {lambdas}. New estimates for lambda will be computed."
                warnings.warn(msg, SyntaxWarning)
                lambdas = (None, None)
        return (self._psi1.estimateThreshold(np.sqrt(overallAccuracy), lambdas[0]),self._psi2.estimateThreshold(np.sqrt(overallAccuracy), lambdas[1]))

    def savePosterior(self, fileName, fileCollisionMethod='rename'):
        """Saves the posterior array over probLambda as a pickle file
        with the specified name.
        Parameters
        ----------
        fileCollisionMethod : string
            Collision method passed to :func:`~psychopy.tools.fileerrortools.handleFileCollision`
        """
        try:
            if os.path.exists(fileName):
                fileName = handleFileCollision(
                    fileName,
                    fileCollisionMethod=fileCollisionMethod
                )
            np.save(fileName, (self._psi1.savePosterior(), self._psi2.savePosterior()))
        except IOError:
            warnings.warn("An error occurred while trying to save the "
                          "posterior array. Continuing without saving...")


# =============================================================================
# Run Adaptive Block
# =============================================================================

def RunAdaptiveBlock (trialFun, numTrials, Xrange, Yrange, Xsteps=100, Ysteps=100, lapseRate=0.0, overallAcc=0.75, logfile=None, info=None, additionalArgs=None):
    """
    Creates an instance of AGRTHandler, then runs the adaptive block.
    User must supply a function to run the individual trial:
        Input: stimulus (tuple of intensities), optional arguments
        Output: tuple of the result (tuple of the response levels) and additional results (tuple or None)
    Maybe optional functions to be called before/after the trial ??
    """
    
    if np.sqrt(overallAcc) >= 1.0 - lapseRate / 2:
        raise RuntimeError('Requested accuracy is too high.')
    elif np.sqrt(overallAcc) <= lapseRate / 2:
        raise RuntimeError('Requested accuracy is too low.')
    
    trialNum = 0
    agrt = AGRTHandler(nTrials=numTrials, lapse=lapseRate, dim1range=Xrange, dim2range=Yrange, dim1steps=Xsteps, dim2steps=Ysteps)
    for intens in agrt:
        trialNum += 1
        result = trialFun(intens, logfile=logfile, handler=agrt, trial=trialNum, info=info, additionalArgs = ['adapt'] if additionalArgs is None else additionalArgs + ['adapt'] if isinstance(additionalArgs, list) else [additionalArgs,'adapt'])
        if result is None: # or whatever value I want to use as the graceful escape code
            raise RuntimeError('Adaptive trial result is None!')
        assert len(result) == 2 and (result[0] == 0 or result[0] == 1) and (result[1] == 0 or result[1] == 1)
        agrt.addResponse(result)
    return agrt.estimateGRTintensities(overallAcc)


# =============================================================================
# GRT Handler
# =============================================================================

class GRTHandler (TrialHandler):
    """
    Same idea as AGRTHandler but for the basic GRT experiment
    May be redundant with standard TrialHandler... double check source code
    Basically implement this as a special case of TrialHandler 
    """

    def __init__(self, nTrials, blockingFactor=2):
        """
        nTrials = total number of trials (4 * number_of_trials_per_stimulus)
        blockFactor = how many times each stimulus should appear in a randomized block (e.g., 2 = shuffle every 8 stimuli (2 of each stimulus))
        """
        
        self.nTrials = nTrials
        self.blockingFactor = blockingFactor
        TrialHandler.__init__(self, trialList = ["(0,0)","(1,0)","(0,1)","(1,1)"] * blockingFactor, 
                              nReps = nTrials / 4 / blockingFactor, 
                              method='random')


# =============================================================================
# Run GRT Block
# =============================================================================

def RunGRTBlock (trialFun, numTrials, stimuli, blockFactor=2, logfile=None, info=None, additionalArgs=None):
    """
    Creates an instance of GRTHandler, then runs the GRT block
    User must supply a function to run the individual trial:
        Input: stimulus, ...
        Output: 
    Maybe optional functions to be called before/after the trial ??
    """
    
    trialNum = 0
    grt = GRTHandler(nTrials=numTrials, blockingFactor=blockFactor)
    
    for mystim in grt:
        stim = eval(mystim) # TrialHandler doesn't like using a tuple of tuples for its trialList so I had to convert them to strings. Convert them back to tuples here.
        trialNum += 1
        result = trialFun((stimuli[0][stim[0]], stimuli[1][stim[1]]), logfile=logfile, handler=grt, trial=trialNum, info=info, additionalArgs = ['main'] if additionalArgs is None else additionalArgs + ['main'] if isinstance(additionalArgs, list) else [additionalArgs,'main'])

        if result is None: # or whatever value I want to use as the graceful escape code
            raise RuntimeError('GRT trial result is None!')
        assert len(result) == 2 and (result[0] == 0 or result[0] == 1) and (result[1] == 0 or result[1] == 1)


# =============================================================================
# Run Adaptive GRT Experiment
# =============================================================================

def RunAdaptiveGRTExperiment (trialFunction, nAdaptiveTrials, nGRTtrials, 
                              dim1range, dim2range, dim1steps=100, dim2steps=100, 
                              lapse=0.0, overallAccuracy=0.75, 
                              adaptiveFunArgs=None, grtFunArgs=None, 
                              blockingFactor=2, logfile=None, info=None):
    """
    Runs the full adaptive + experiment using the respective handlers
    """
    
    RunGRTBlock(trialFun=trialFunction, numTrials=nGRTtrials, 
                stimuli=RunAdaptiveBlock(trialFun=trialFunction, numTrials=nAdaptiveTrials, 
                                         Xrange=dim1range, Yrange=dim2range, Xsteps=dim1steps, 
                                         Ysteps=dim2steps, lapseRate=lapse, overallAcc=overallAccuracy, 
                                         additionalArgs=adaptiveFunArgs, logfile=logfile, info=info), 
                blockFactor=blockingFactor, additionalArgs=grtFunArgs, logfile=logfile, info=info)


# =============================================================================
# Subject Model
# =============================================================================

def GRTSubjectModel (stimulus, logfile=None, handler=None, trial=None, info=None, additionalArgs=None):
    """
    Given stimulus coordinates and model parameters, returns the response probabilities for stimulus index R.
    If R is not supplied, all response probabilities are returned.

    x := value on stimulus dimension x
    y := value on stimulus dimension y
    m := 2 element np.array specifying the mean spreading for each dimension
    n := 2 element np.array specifying the var spreading for each dimension
    rho := 4 element np.array specifying the correlation for each quadrant
    sig := 2 element np.array specifying the standard deviation for each dimension
    delta := 2 element np.array specifying the decision bound intercept for each dimension
    epsilon := 2 element np.array specifying the decision bound slopes for each dimension
    """
    
    x, y = stimulus
    delta, epsilon, sig, m, n, rho, lapse = additionalArgs[0] # R=[0,1,2,3]
    
    # calculate decision bounds based on where we are in stimulus space
    tau = delta + epsilon * np.array([y,x])
    # create indicator for which quadrant we're in
    zeta = np.array([1 if x > tau[0] else -1, 1 if y > tau[1] else -1])
    # calculate mean vector
    mymean = np.array([x,y]) + np.flip(zeta) * m * np.array([x,y])
    # calculate sd vector
    mysd = sig + np.flip(zeta) * n * sig
    # pick the rho that corresponds to our quadrant
    if sum(zeta) == 2:
        myrho = rho[3]
    elif sum(zeta) == -2:
        myrho = rho[0]
    elif zeta[0] > 0:
        myrho = rho[2]
    else:
        myrho = rho[1]
    # create the cov matrix
    mycov = np.array([[mysd[0]**2, myrho * mysd[0] * mysd[1]], [myrho * mysd[0] * mysd[1], mysd[1]**2]])
    
    if random.random() < lapse:
        result = random.choices(((0,0),(0,1),(1,0),(1,1)))
    else:
        result = random.choices(((0,0),(0,1),(1,0),(1,1)), 
                                np.array([mvnun(np.array([-np.inf, -np.inf]), tau, mymean, mycov)[0], 
                                          mvnun(np.array([-np.inf, tau[1]]), np.array([tau[0], np.inf]), mymean, mycov)[0],
                                          mvnun(np.array([tau[0], -np.inf]), np.array([np.inf, tau[1]]), mymean, mycov)[0],
                                          mvnun(tau, np.array([np.inf, np.inf]), mymean, mycov)[0]]))
    [result] = result # random.choices returns a list, but since I only want the single response tuple, unpack the returned list
    
    # log results
    if logfile is not None:
        flatten = lambda l: sum(map(flatten,l),[]) if isinstance(l,list) else [l]
        logstr = flatten([i for i in [info, additionalArgs, trial, stimulus] if i is not None])
        logstr = "\t".join(map(str, logstr + [result, 'rt', handler.estimateLambda() if isinstance(handler, AGRTHandler) else 'NA'])) + "\n"
        logfile.write(logstr)
    
    return result

