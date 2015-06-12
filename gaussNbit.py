import numpy as np
import scipy
import scipy.special
import scipy.stats
import math

#-------- erf(x/sqrt(2))
def erf2(x):
    return scipy.special.erf(x / math.sqrt(2.0))

#-------- probabilities among uniformly-spaced threshold voltages
def probNbit(param, levelNum):
    #    param: 2-element vector, param[0] is the voltage interval (V0/sigma) and param[1] is bias (mu/sigma), scaled by SD
    volt = param[0]* np.linspace( -levelNum/2+1, levelNum/2-1, levelNum - 1) - param[1]
    prob = np.r_[1.0 + erf2(volt[0]), erf2(volt[1:(levelNum-1)])-erf2(volt[0:(levelNum-2)]), 1.0 - erf2(volt[levelNum - 2])]
    return 0.5* prob

#-------- probabilities among non-uniformly-spaced given threshold voltages
def probNbitThresh(param, levelNum, thresh):
    #    param: 2-element vector, param[0] is the voltage interval (V0/sigma) and param[1] is bias (mu/sigma), scaled by SD
    volt = param[0]* thresh - param[1]
    prob = np.r_[1.0 + erf2(volt[0]), erf2(volt[1:(levelNum-1)])-erf2(volt[0:(levelNum-2)]), 1.0 - erf2(volt[levelNum - 2])]
    return 0.5* prob

#-------- Initial Parameters of Gaussian Fit
def initGaussNbit( prob, levelNum ):
    Vweight = np.linspace(-levelNum/2 + 0.5, levelNum/2 - 0.5, levelNum)    # Voltage for each level
    Pweight = Vweight**2                                                    # Power of each level
    Average = np.dot(Vweight, prob)
    SD        = math.sqrt(np.dot(Pweight, prob) - Average**2)
    return 1.0/SD, Average/SD
#
#-------- Gaussian parameters derermined by bit distribution
# gaussNbit(): a function to estimate standard deviation and bias voltage using the level histogram in quantized voltages.
# It assumes that the threshold level of the digitizer is spaced uniformly, and the input signal follows Gaussian distribution.
# The output values are:
#    standard deviation, scaled by the threshold voltage
#    bias voltage, scaled by the threshold voltage
#    stnadard error of the standard deviation
#    stnadard error of the bias voltage
#
# Note that the (analog) power will be proportional to 1/gaussNbit()[0][0]**2
def gaussNbit( Nsample, levelNum ):
    #    Nsample:    Histogram of each quantization level, as an array with [levelNum] elements
    #    levelNum:    Number of quantization levels (for ALMA, levelNum = 8)
    #
    if(len(Nsample) < levelNum):
        return    -1
    #
    #-------- Initial Statistics
    Nsample = np.array(Nsample, dtype=float)
    totalSample = np.sum(Nsample[0:levelNum])    # Total number of samples
    prob = Nsample[0:levelNum]/totalSample        # Measured probability in each state
    weight = Nsample[0:levelNum] / (1.0 - prob)**2    # Weight following biniminal distribution
    W = np.diag( weight )                            # Weight matrix
    SqPi2 = 1.0 / math.sqrt(2.0* math.pi)
    a = initGaussNbit(prob, levelNum)            # Initial parameters
    #
    #-------- Loop for iteration
    niter = 0
    while(niter < 10):
        resid = prob - probNbit( a, levelNum )    # Residual of measured histogram from expected
        #
        #-------- Factors used in the partial matrix
        erfDeriv = np.r_[0.0, np.exp(-0.5* (np.linspace(-levelNum/2+1, levelNum/2-1, levelNum-1)* a[0] - a[1])**2), 0.0]
        WerfDeriv = np.linspace(-levelNum/2, levelNum/2, levelNum+1)* erfDeriv
        #
        #-------- Partial Matrix
        p = np.array([SqPi2* np.diff(WerfDeriv),  -SqPi2* np.diff(erfDeriv)])
        t_p_W = np.dot(p, W)
        #-------- Solution
        solution = np.linalg.inv(np.dot(t_p_W, p.T))
        correction = np.dot(solution, np.dot(t_p_W, resid))
        #-------- Correction
        a = a + correction
        #-------- Converged?
        if( np.dot(correction, correction) < 1.0e-20 ): break
        niter = niter + 1
    #
    return a, np.diag(solution)
#
#-------- Gaussian parameters estimation with non-uniformly spaced given thresholds
def gaussNbitThresh( Nsample, levelNum, thresh ):
    #    Nsample:    Histogram of each quantization level, as an array with [levelNum] elements
    #    levelNum:    Number of quantization levels (for ALMA, levelNum = 8)
    #   thresh:     Given threshold voltages
    #
    if(len(Nsample) < levelNum):
        return    -1
    #
    #-------- Initial Statistics
    Nsample = np.array(Nsample, dtype=float)
    totalSample = np.sum(Nsample[0:levelNum])    # Total number of samples
    prob = Nsample[0:levelNum]/totalSample        # Measured probability in each state
    weight = Nsample[0:levelNum] / (1.0 - prob)**2    # Weight following biniminal distribution
    W = np.diag( weight )                            # Weight matrix
    SqPi2 = 1.0 / math.sqrt(2.0* math.pi)
    a = initGaussNbit(prob, levelNum)            # Initial parameters
    #
    #-------- Loop for iteration
    niter = 0
    while(niter < 10):
        resid = prob - probNbitThresh( a, levelNum, thresh )    # Residual of measured histogram from expected
        #
        #-------- Factors used in the partial matrix
        erfDeriv = np.r_[0.0, np.exp(-0.5* (thresh* a[0] - a[1])**2), 0.0]
        WerfDeriv = np.r_[0.0, thresh, 0.0]* erfDeriv
        #
        #-------- Partial Matrix
        p = np.array([SqPi2* np.diff(WerfDeriv),  -SqPi2* np.diff(erfDeriv)])
        t_p_W = np.dot(p, W)
        #-------- Solution
        solution = np.linalg.inv(np.dot(t_p_W, p.T))
        correction = np.dot(solution, np.dot(t_p_W, resid))
        #-------- Correction
        a = a + correction
        #-------- Converged?
        if( np.dot(correction, correction) < 1.0e-20 ): break
        niter = niter + 1
    #
    return a, np.diag(solution)
#
#-------- Estimation of threshold voltages 
def threshLevel( Nsample ):
    levelNum = len(Nsample)
    gaussParam = gaussNbit( Nsample, levelNum)[0]
    Nsample = np.array(Nsample, dtype=float); totalSample = np.sum(Nsample)
    probReal = Nsample[0:levelNum]/totalSample
    cumReal = probReal[0:(levelNum-1)]
    for level_index in range(levelNum - 1):
        cumReal = cumReal + np.r_[ np.zeros(level_index+1), probReal[0:-level_index-2]]
    #
    thresh = (scipy.stats.norm.ppf(cumReal) + gaussParam[1]) / gaussParam[0]
    return thresh
#
#-------- Estimation of threshold voltages 
def threshLevelRef( Nsample, refThresh ):
    levelNum = len(Nsample)
    gaussParam = gaussNbitThresh( Nsample, levelNum, refThresh)[0]
    Nsample = np.array(Nsample, dtype=float); totalSample = np.sum(Nsample)
    probReal = Nsample[0:levelNum]/totalSample
    cumReal = probReal[0:(levelNum-1)]
    for level_index in range(levelNum - 1):
        cumReal = cumReal + np.r_[ np.zeros(level_index+1), probReal[0:-level_index-2]]
    #
    thresh = (scipy.stats.norm.ppf(cumReal) + gaussParam[1]) / gaussParam[0]
    return thresh
#
#-------- Power derived from uniform thresholds
def uniPower( sigma ):
    weight = arange(-3.5, 4.5, 1)
    return np.sum( probNbit( np.array([sigma, 0.0]), 8)* weight**2)
#
#-------- Power derived from uniform thresholds
def threshPower( sigma, thresh ):
    weight = arange(-3.5, 4.5, 1)
    return np.sum( probNbitThresh( np.array([sigma, 0.0]), 8, thresh)* weight**2)
#
