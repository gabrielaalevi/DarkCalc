#This module has auxiliary functions
from scipy.special import kn, zeta
from numpy import pi,sqrt,exp
m_planck = 1.2209 * 10**(19) #reduced planck mass, in GeV

def gstar(T: float) -> float:
    """
    Calculates the number of relativistic degrees of freedom for a given temperature T (in GeV)
    assuming only the SM degrees of freedom.

    :param T: thermal bath temperature (in GeV)

    """
    
    if T >= 34:
        g_star = 106.75
    elif T>= 14:
        g_star = 90
    elif T>= 4.43:
        g_star = 85
    elif T>= 1.75:
        g_star = 80
    elif T >= 0.7:
        g_star = 70
    elif T>= 0.28:
        g_star = 60
    elif T>= 0.2:
        g_star = 50
    elif T >= 0.1:
        g_star = 18
    elif T >= 0.0132:
        g_star = 10
    else:
        g_star = 3
    
    return g_star


def H(T: float) -> float:
    """
    Calculates the Hubble rate for a given temperature T (in GeV) assuming a thermal bath
    dominated by radiation and SM degrees of freedom.

    :param T: thermal bath temperature (in GeV)
    """
    
    g_star = gstar(T)
    h = 2*sqrt(g_star*pi**3/45)*T**2/m_planck
    return h


def S(T: float) -> float:
    """
    Calculates the entropy density s for a given temperature T (in GeV)

    :param T: thermal bath temperature (in GeV)
    """
    
    g_star = gstar(T) #g_star is the number of relativistic degrees of freedom
    s = (2*pi**2/45)*g_star*T**3

    return s


def Neq(T: float, m: float, g: int) -> float:
    """
    Calculates the equilibrium number density at a temperature T (in GeV) for a certain particle, 
    using Fermi-Dirac or Bose-Einstein Statistics.

    :param T: thermal bath temperature (in GeV)
    :param m: Particle mass
    :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

    """
    #m represents particle's mass, T is temperature in GeV, g represents the particle's degrees of freedom, 
    #and s represents entropy density.
    

    xm = m/T     
    if xm > 10: #non-relativistic regime
        coeffs = [1,15./8.,105./128.,-315./1024.,10395./32768.]
        neq = (m**3*exp(-xm)/(xm**(3/2)))*(1/(2*sqrt(2)*pi**(3/2)))
        neq = neq*sum([c/xm**i for i,c in enumerate(coeffs)])
    elif xm > (2./3.): #semi-relativistic regime
        neq = m**3*(1/(2*pi**2))*kn(2, xm)/xm
    else: #relativistic regime
        neq = (zeta(3)/pi**2)*T**3 # bosons
        if g < 0:
           neq = (3/4)*neq # fermions
            
    neq = abs(g)*neq
    
    return neq

def Yeq(T: float, m: float, g: int) -> float:
    """
    Calculates the equilibrium yield at a temperature T (in GeV) for a certain particle, 
    using Fermi-Dirac or Bose-Einstein Statistics.

    :param T: thermal bath temperature (in GeV)
    :param m: Particle mass
    :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

    """
    
    
    s = S(T) #entropy density
    neq = Neq(T,m,g)
    yeq = neq/s
    return yeq


def dSdT(T: float) -> float: 
    """
    Calculates the variation of entropy with respect to T.

    :param T: thermal bath temperature (in GeV)
    """
    
    g_star = gstar(T) #g_* is the number of relativistic degrees of freedom    
    dsdT = (6*pi**2/45)*g_star*T**2
    
    return dsdT

def dSdx(x: float, mDM: float) -> float:
    """
    Calculates the variation of entropy with respect to x=mDM/T.

    :param T: thermal bath temperature (in GeV)
    """
    
    T = mDM/x    
    g_star = gstar(T) #g_* is the number of relativistic degrees of freedom    
    dsdT = (6*pi**2/45)*g_star*T**2
    dsdx = (-mDM/x**2)*dsdT
    
    return dsdx
