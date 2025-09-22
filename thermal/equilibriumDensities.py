#This module has auxiliary functions
from scipy.special import kn, zeta
from numpy import pi,sqrt,exp
import numpy as np
from tools.logger import logger
m_planck = 1.2209 * 10**(19) #reduced planck mass, in GeV

def gstar(T: float) -> float:
    """
    Calculates the number of relativistic degrees of freedom for a given temperature T (in GeV)
    assuming only the SM degrees of freedom.

    :param T: thermal bath temperature (in GeV)

    """
    
    if T >= 59:
        g_star = 106.75
    elif T>= 49:
        g_star = 103
    elif T >= 40:
        g_star = 102
    elif T>= 34:
        g_star = 100
    elif T>= 23:
        g_star = 95
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


def Neq_n(n: int, T: float, m: float, g: int) -> float:
    """
    Calculates the equilibrium number density at a temperature T (in GeV) for a certain particle, 
    using Fermi-Dirac or Bose-Einstein Statistics.

    :param T: thermal bath temperature (in GeV)
    :param m: Particle mass
    :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

    """

    xm = m/T
    fac = abs(g)
    if g < 0: 
        fac = fac*(-1)**n # factor for fermions
    
    if xm > 10: #non-relativistic regime
        coeffs = np.array([1,15./8.,105./128.,-315./1024.,10395./32768.])
        fac = fac*(m**3*exp(-xm*(1+n))/(xm**(3/2)))*(1/(2*sqrt(2)*pi**(3/2)))
        terms = np.array([(xm**i)*(1.0+n)**(3./2.+i) for i in range(len(coeffs))])
        neq_n = fac*sum(coeffs/terms)
    elif xm > 0.15: #semi-relativistic regime
        neq_n = fac*m**3*(1/(2*(1+n)*pi**2))*kn(2, xm*(1+n))/xm
    else: #relativistic regime
        if n == 0:
            neq_n = abs(g)*(zeta(3)/pi**2)*T**3 # bosons
            if g < 0:
                neq_n = (3/4)*neq_n # fermions
        else:
            neq_n = 0.0
    
    return neq_n

def Neq(T: float, m: float, g: int, nmax : int = 5) -> float:
    """
    Calculates the equilibrium number density at a temperature T (in GeV) for a certain particle, 
    using Fermi-Dirac or Bose-Einstein Statistics.

    :param T: thermal bath temperature (in GeV)
    :param m: Particle mass
    :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions
    :param nmax: Maximum order for corrections to the Maxwerll-Boltzmann distribution 
                 (nmax = 0 is equivalent to assuming Maxwell-Boltzmann for the non-relativistic regime)

    """

    neq = 0.0
    for n in range(0,nmax+1):
        neq += Neq_n(n,T,m,g)
    
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

def Req(T: float, m1: float, g1: int, m2 : float, g2: int) -> float:
    """
    Calculates the ratio of equilibrium number densities at a temperature T (in GeV), i.e. Neq(1)/Neq(2) 
    using Maxwell-Boltzmann statitics.

    :param T: thermal bath temperature (in GeV)
    :param m: Particle mass
    :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

    """

    xm1 = m1/T
    xm2 = m2/T
    # If there is no strong Boltzmann suppression, directly compute the ratio
    if max(xm1,xm2) < 10.0:
        return Neq(T,m1,g1)/Neq(T,m2,g2)

    coeffs = np.array([1,15./8.,105./128.,-315./1024.,10395./32768.])
    if xm1 == xm2:
        mdiff_fac = 1.0
    else:
        dm = (xm1-xm2)
        if xm1 == 0.0:
            rm = 1/xm2**(3/2)
        elif xm2 == 0.0:
            rm = xm1**(3/2)
        else:
            rm = (xm1/xm2)**(3/2)
        mdiff_fac = rm*np.exp(-dm)

    xi = []
    for xm,g in zip([xm1,xm2],[g1,g2]):
        if xm > 10: #non-relativistic regime         
            fac = (1/(2*sqrt(2)*pi**(3/2)))       
            terms = np.array([xm**i for i in range(len(coeffs))])
            xi.append(fac*sum(coeffs/terms))
        elif xm > 0.15: #semi-relativistic regime
            fac = 1/(2*pi**2)
            xi.append(fac*np.sqrt(xm)*np.exp(xm)*kn(2, xm))
        else: #relativistic regime
            if g > 0:
                xi.append(zeta(3)/pi**2)
            else:
                xi.append((3/4)*zeta(3)/pi**2)
    
    if xi[1] <= 0.0:
        req = np.nan
    elif xi[0] <= 0.0:
        req = 0.0
    else:
        req = (xi[0]/xi[1])*mdiff_fac*abs(g1)/abs(g2)

    if np.isinf(req) or np.isnan(req):
        logger.error(f"Trying to compute ratios of number densities at T = {T:1.2e} lead to infinities. Try to increase the final temperature.")
        raise ValueError()
    
    return req


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
