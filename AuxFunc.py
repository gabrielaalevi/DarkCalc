#This module has auxiliary functions
import math
from scipy.special import kn, zeta
m_planck = 2.4 * 10**(18) #reduced planck mass, in GeV

#calculates the number of relativistic degrees of freedom for temperature T (in GeV), using solely the SM degrees of freedom
def degrees_of_freedom(T):
    if T >= 0.3:
        g_star = 100
    elif T >= 0.1:
        g_star = 10
    elif 0.1 > (T):
        g_star = 3
    return g_star

#calculates the hubble rate H based on the temperature T (in GeV)
def hubblerate(T):
    g_star = degrees_of_freedom(T) #g_star is the number of relativistic degrees of freedom
    hubble = ((math.pi**2)/90)**(1/2) * (g_star **(1/2)) * T**2 / m_planck #hubble constant, varying with temperature
    return hubble

#calculates the entropy density s, depending on the temperature T (in GeV)
def entropydensity(T):
    g_star = degrees_of_freedom(T) #g_star is the number of relativistic degrees of freedom
    s = (2/45) * (math.pi**2) * (g_star) * T**3 #entropy density
    return s

#calculates the equilibrium yield for a certain particle, using Fermi-Dirac or Bose-Einstein Statistics.
def equilibrium_yield(m, T, g):
    #m represents particle's mass, T is temperature in GeV, g represents the particle's degrees of freedom, 
    #and s represents entropy density.
    g_star = degrees_of_freedom(T)#g_* is the number of relativistic degrees of freedom
    s = (2/45) * (math.pi**2) * (g_star) * T**3 #entropy density
    if (m/T)>10: #non-relativistic regime
        Y_eq = (g * (m * T/ (2 * math.pi))**(3/2) * math.exp(-m/T))/s
    elif (m/T)>(2/3): #semi-relativistic regime
        Y_eq = (45 * g * (m/T)**2 * kn(2,(m/T)))/(4 * math.pi**4 * g_star)
    elif (m/T)<=(2/3): #relativistic regime
        if g > 0:
            Y_eq = ((3/4) * (zeta(3)/math.pi**2) * g * (T)**3)/s
        if g <= 0:
            Y_eq = ((zeta(3)/math.pi**2) * g * (T)**3)/s
    return Y_eq

#calculates the variation of entropy with x, where x = m_DM/T
def dsdx(x,m): #m is the mass of the DM particle
    g_star = degrees_of_freedom(m/x) #g_* is the number of relativistic degrees of freedom
    ds_dx = (-6 * math.pi**2 * g_star * m**3)/(45 * x**4)
    return ds_dx