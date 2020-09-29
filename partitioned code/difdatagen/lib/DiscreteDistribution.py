import random
import math
 
class DiscreteDistribution():
    """ Discrete probability distribution generating class
    """
    def __init__(self):
       return

    @classmethod
    def getExpoPdf(self,lam,distlen,add=1):
       pdf = {}
       mylambda = lam[0]
       for index in xrange(distlen):
           val = math.e**(-1.0*mylambda*index) - math.e**(-1.0*mylambda*(index+1))
           pdf[index + add] = val
       return pdf
    
    @classmethod
    def getWeibullPdf(self,params,distlen,add=1):
       pdf = {} 
       scale,shape = params
       for index in xrange(distlen):
           upperval = 1.0-math.e**(-1.0*((float(index+1)/scale)**shape))
           lowerval = 1.0-math.e**(-1.0*((float(index)/scale)**shape))
           pdf[index+add] = upperval-lowerval
       return pdf

    @classmethod
    def getRayleighPdf(self,paramscale,distlen,add=1):
       paramscale = paramscale[0]
       pdf = {} 
       for index in xrange(distlen):
           upperval = 1.0-math.e**(-1.0*(((float(index+1)**2)/(2.0*(paramscale**2)))))
           lowerval = 1.0-math.e**(-1.0*(((float(index)**2)/(2.0*(paramscale**2)))))
           pdf[index+add] = upperval-lowerval
       return pdf

    @classmethod
    def getLognormalPdf(self,params,distlen,add=1):
       pdf = {} 
       mu, sigma = params
       for index in xrange(distlen):
           useindex = index + 0.000000001
           upperval = 0.5 + (0.5*math.erf(float(math.log(useindex+1)-mu)/math.sqrt(2.0*(sigma**2))))
           lowerval = 0.5 + (0.5*math.erf(float(math.log(useindex)-mu)/math.sqrt(2.0*(sigma**2))))
           pdf[index+add] = upperval - lowerval
       return pdf

    @classmethod 
    def genPdf(self,partdict):
        """ generates pdf of given part distribution
        """
        pdf,cur = {},1.0
        for key in sorted(partdict.keys()):
            pdf[key] = cur * partdict[key]
            cur = pdf[key]
        return pdf     
            
    @classmethod
    def genPartDist(self,dist,distparam,mode,sprob=1.0,distlen=20):
       """generates discrete partition distribution
       Args:
          dist: distribution
          distparam: distribution parameters
          mode: dist mode
          distlen: distribution length
          sprob: spreading probability
       Returns:
          time2ratio: ratio per time of given dist
       """
       assert mode in ["normal", "reverseCdf", "normalCdf"]
       assert dist in ["expo","rayleigh","weibull","lognormal"]
       pdffunc = "get{0}Pdf".format(dist.capitalize())
       method = getattr(self, pdffunc)
       time2prob = method(distparam,distlen)
       assert sum(time2prob.values()) <= 1.000000000001
       maxtime = max(time2prob.keys())
       time2prob = {time: prob * sprob for (time,prob) in time2prob.items()}
       if mode == "normal":
          rightsize = dict(time2prob) 
       elif mode == "reverseCdf":
          rightsize = {}
          mysum = 0.0
          for index in xrange(1,maxtime+1):
              mysum += time2prob[index]
              rightsize[index] = 1.0 - mysum
              if rightsize[index] < 0:
                 rightsize[index] = 0.0
              assert rightsize[index] >= 0
       elif mode == "normalCdf":
          rightsize = {}
          mysum = 0.0
          for index in xrange(1,maxtime+1):
              mysum += time2prob[index]
              rightsize[index] = mysum
       time2ratio={}   
       time2ratio[1] = rightsize[1] 
       for index in xrange(2,distlen + 1):
           if rightsize[index-1] == 0:
              ratio = 0.0
           else:   
              ratio = float(rightsize[index]) / rightsize[index-1] 
           time2ratio[index] = ratio
       return time2ratio
