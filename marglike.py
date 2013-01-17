#!/cygdrive/f/Python27/python.exe

import sys
import csv
import numpy
from scipy.stats.mstats import mquantiles
from scipy.optimize import curve_fit

class MarginalLikelihood(object):
  def __init__(self,
    filename,
    burnin=None,
    delimiter='\t',
    chainIdColumn='chainId',
    logLikelihoodColumn='logLikelihood',
    heatColumn='heat',
    heatPower=1.0
  ):
    # Read data into a lists
    chainIds = list()
    heats = list()
    logLikelihoods = list()
    
    f = open(filename)
    reader = csv.DictReader(f, delimiter=delimiter)
    for row in reader:
      chainIds.append(int(row[chainIdColumn]))
      if 'heat' in row:
        heats.append(float(row[heatColumn]))
      logLikelihoods.append(float(row[logLikelihoodColumn]))
    f.close()
    
    # Get number of chains and check length of log-likelihoods
    self.chainCount = max(chainIds) + 1
    if not (len(logLikelihoods) % self.chainCount == 0):
      raise Exception('%d % %d'.format(len(logLikelihoods), self.chainCount),
        'Number of log-likelihoods is not a multiple of the number of chains'
      )
    
    self.heats = numpy.zeros(self.chainCount, dtype=float)
    self.heats.fill(-1.0)
    heatCount = 0
    if len(heats) == len(chainIds):
      self.heatPower = None
      for heat, chainId in zip(heats, chainIds):
        if heats[chainId] == -1.0:
          self.heats[chainId] = heat
          heatCount += 1
        if heatCount == len(chainIds):
          break
    else:
      self.heatPower = heatPower
      for i in range(self.chainCount):
        self.heats[i] = (float(self.chainCount - 1 - i) /
          (self.chainCount - 1))**self.heatPower
    
    assert self.heats[0] == 1.0
    assert self.heats[self.chainCount - 1] == 0.0
    
    # Construct data matrix
    self.data = numpy.zeros((self.chainCount, len(logLikelihoods) / self.chainCount))
    for i, (chainId, logLikelihood) in enumerate(zip(chainIds, logLikelihoods)):
      self.data[chainId, i/self.chainCount] = logLikelihood
    
    self.burnin = burnin if burnin != None else self.data.shape[1] / 2
  
  def bootstrap(self, func, sampleCount, quantiles=(0.025, 0.975)):
    estimates = numpy.zeros(sampleCount)
    for i in range(sampleCount):
      sample = numpy.zeros((self.chainCount, self.data.shape[1] - self.burnin))
      for chainId in range(self.chainCount):
        size = self.data.shape[1] - self.burnin
        sample[chainId,:] = numpy.random.choice(
          self.data[chainId,self.burnin:], size=size
        )
      estimates[i] = func(data=sample)
    quantiles = mquantiles(estimates, quantiles)
    
    return estimates, quantiles  
    
  def trapezoidRuleEstimate(self, data=None):
    if data == None:
      data = self.data[:,self.burnin:]
    
    means = numpy.mean(data, axis=1)
    ml = 0.0
    for i in range(self.chainCount - 1):
      ml += (self.heats[i] - self.heats[i+1]) * (means[i+1] + means[i]) / 2.0
    return ml
  
  def trapezoidRuleBootstrap(self, sampleCount, quantiles=(0.025, 0.975)):
    return self.bootstrap(self.trapezoidRuleEstimate, sampleCount, quantiles)
  
  def curveFitEstimate(self, data=None):
    def curve(t, l1, l0, p):
      return l1 - (l1 - l0) * (1.0 - t)**p
    
    if data == None:
      data = self.data[:,self.burnin:]
    means = numpy.mean(data, axis=1)
    
    params, cov = curve_fit(curve, self.heats, means, p0=(means[0], means[-1], 1))
    
    l1 = params[0]
    l0 = params[1]
    p = params[2]
    
    ml = l1 - (l1 - l0) / (p + 1)
    return ml
  
  def curveFitBootstrap(self, sampleCount, quantiles=(0.025, 0.975)):
    return self.bootstrap(self.curveFitEstimate, sampleCount, quantiles)
  

def runMargLike(filename, heatPower):
  ml = MarginalLikelihood(filename, heatPower=heatPower)
  print 'trapezoid-rule estimate: ', ml.trapezoidRuleEstimate()
  estimates, quantiles = ml.trapezoidRuleBootstrap(1000)
  print 'trapezoid-rule bootstrap 95% interval: ', quantiles
  
  print 'curve-fit estimate: ', ml.curveFitEstimate()
  estimates, quantiles = ml.curveFitBootstrap(1000)
  print 'curve-fit bootstrap 95% interval: ', quantiles

if __name__ == '__main__':
    runMargLike(sys.argv[1], float(sys.argv[2]))
