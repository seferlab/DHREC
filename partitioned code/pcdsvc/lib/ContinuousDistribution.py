import random
import math
import scipy.special

class ContinuousDistribution():
    """ Continuous probability distribution generating class
    """
    def __init__(self):
        return

    @classmethod
    def getDistMean(self,dist,param):
        """ returns distribution mean
        """
        if dist == "expo":
           mean  = 1.0 / param[0]
        elif dist == "rayleigh":
           mean  = param[0] * math.sqrt(math.pi/2.0)
        elif dist == "weibull":
           lam, k = float(param[0]), param[1]
           mean = lam * scipy.special.gamma(1.0 + (1.0 / k))
        return mean
    
    @classmethod
    def getExpoPdf(self,lam,x,sprob):
        lam = lam[0]
        return sprob*lam*math.exp(-1.0*lam*x)
    @classmethod
    def getPowerlawPdf(self,alpha,x,sprob):
        alpha = alpha[0]
        return sprob*alpha*math.pow(x,(-1.0*alpha)-1.0)
    @classmethod
    def getLognormalPdf(self,params,x,sprob):
        mu, sigma = params
        return sprob*(1.0/(math.sqrt(2.0*math.pi)*x*sigma))*math.exp(-0.5*((math.log(x)-mu)**2)/(sigma**2))
    @classmethod 
    def getWeibullPdf(self,params,x,sprob):
        lam, k = float(params[0]), params[1]
        print "info"
        print lam
        print k
        return sprob*(k/lam)*((x/lam)**(k-1))*math.exp(-1.0*((x/lam)**k))
    @classmethod 
    def getRayleighPdf(self,sigma,x,sprob):
        sigma = sigma[0]
        return sprob*(x/(sigma**2))*math.exp((-0.5*x*x)/(sigma**2))
    @classmethod 
    def getExpoCdf(self,lam,x,sprob):
        lam = lam[0]
        return sprob*(1.0-math.exp(-1.0*lam*x))
    @classmethod 
    def getPowerlawCdf(self,alpha,x,sprob):
        alpha = alpha[0]
        return sprob*(1.0-math.pow(x,-1.0*alpha))
    @classmethod
    def getLognormalCdf(self,params,x,sprob):
        mu,sigma = params
        return sprob*0.5*(1.0+math.erf((math.log(x)-mu)/(sigma*math.sqrt(2))))
    @classmethod
    def getWeibullCdf(self,params,x,sprob):
        lam,k = float(params[0]), params[1]
        return sprob*(1.0-math.exp(-1.0*((x/lam)**k)))
    @classmethod 
    def getRayleighCdf(self,sigma,x,sprob):
        sigma = sigma[0]
        return sprob*(1.0-math.exp((-0.5*x*x)/(sigma**2)))
    @classmethod
    def getExpoRandNum(self,lam):
        lam = lam[0]
        return math.log(random.random())/(-1.0*lam)
    @classmethod
    def getPowerlawRandNum(self,alpha):
        alpha = alpha[0]
        return math.pow(random.random(),-1.0/alpha)
    @classmethod
    def getRayleighRandNum(self,sigma):
        sigma = sigma[0]
        num = math.log(random.random())*-2.0*(sigma**2)
        return math.sqrt(num)
    @classmethod
    def getWeibullRandNum(self,params):
        lam,k = float(params[0]), params[1] 
        num = math.log(random.random())*-1.0
        return math.pow(num,1.0/k)*lam
    @classmethod
    def getLognormalRandNum(self,params):
        mu,sigma = params
        return math.exp(mu+(sigma*random.gauss(0,1.0)))  

    @classmethod
    def genContRandNum(self,dist,params,sprob=1.01):
        """generates continous random variable
        Args:
          dist: distribution 
          params: dist. parameters
          sprob: s
        Returns:
          randnum: returns continuous random number
        """ 
        assert dist in ["expo","powerlaw","rayleigh","weibull","lognormal"] 
        if random.random() >= sprob:
           return None
        pdffunc = "get{0}RandNum".format(dist.capitalize())
        method = getattr(self, pdffunc)
        return method(params)

    @classmethod
    def getContPdf(self,dist,params,x,sprob=1.0):
       """ returns pdf(x)
       Args:
          dist: distribution 
          params: dist. parameters
          x: x point
          sprob: spreading probability
       Returns:
          f(x): returns continous random number
       """
       assert dist in ["expo","powerlaw","rayleigh","weibull","lognormal"]
       pdffunc = "get{0}Pdf".format(dist.capitalize())
       method = getattr(self, pdffunc)
       return method(params,x,sprob)

    @classmethod
    def getContCdf(self,dist,params,x,sprob=1.0):
       """ returns cdf(x)
       Args:
          dist: distribution 
          params: dist. parameters
          x: x point
          sprob: spreading probability
       Returns:
          f(x): returns continous random number
       """
       assert dist in ["expo","powerlaw","rayleigh","weibull","lognormal"]
       pdffunc = "get{0}Cdf".format(dist.capitalize())
       method = getattr(self, pdffunc)
       return method(params,x,sprob)
