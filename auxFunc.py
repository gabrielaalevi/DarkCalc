#This module has auxiliary functions
import math
from scipy.special import kn, zeta
import numpy as np
m_planck = 2.4 * 10**(18) #reduced planck mass, in GeV

#calculates the number of relativistic degrees of freedom for temperature T (in GeV), using solely the SM degrees of freedom
def degrees_of_freedom(x, mDM):
    T = mDM/x
    if T >= 1:
        g_star = 106.75
    elif T >= 0.04:
        g_star = 61
    elif T >= 0.0005:
        g_star = 10
    else:
        g_star = 3
    return g_star

#calculates the hubble rate H based on the temperature T (in GeV)
def hubblerate(x, mDM):
    T = mDM/x
    g_star = degrees_of_freedom(x, mDM) #g_star is the number of relativistic degrees of freedom
    hubble = (((math.pi**3)/45)**(1/2) * (g_star **(1/2)) * 2 * T**2) / m_planck #hubble constant, varying with temperature
    return hubble

#calculates the entropy density s, depending on the temperature T (in GeV)
def entropydensity(x, mDM):
    T = mDM/x
    g_star = degrees_of_freedom(x, mDM) #g_star is the number of relativistic degrees of freedom
    s = (2/45) * (math.pi**2) * (g_star) * T**3 #entropy density
    return s

#calculates the equilibrium yield for a certain particle, using Fermi-Dirac or Bose-Einstein Statistics.
def equilibrium_yield(m, x, mDM, g):
    #m represents particle's mass, T is temperature in GeV, g represents the particle's degrees of freedom, 
    #and s represents entropy density.
    T = mDM/x
    g_star = degrees_of_freedom(x, mDM)#g_* is the number of relativistic degrees of freedom
    s = entropydensity(x, mDM) #entropy density
    if (m/T)>10: #non-relativistic regime
        Y_eq = (45/(4 * np.pi**4)) * (abs(g)/g_star) * (m/T)**(1.5) * np.sqrt(np.pi/2) * np.exp(-m/T) * (1 + (15/(8 * (m/T))) + 105/(128 * (m/T)**2) - (315/(1024 * (m/T)**3) )) #(abs(g) * (m * T/ (2 * math.pi))**(3/2) * math.exp(-m/T))/s
    if (m/T)>(2/3): #semi-relativistic regime
        Y_eq = ((abs(g)/(2 * np.pi**2)) * ( T * m**2)* kn(2, m/T)/s) #(45 * abs(g) * (m/T)**2 * kn(2,(m/T)))/(4 * math.pi**4 * g_star)
    elif (m/T)<=(2/3): #relativistic regime
        if g > 0: #fermions
           Y_eq = ((3/4) * (zeta(3)/math.pi**2) * abs(g) * (T)**3)/s
        if g < 0: #bosons
            Y_eq = ((zeta(3)/math.pi**2) * abs(g) * (T)**3)/s
    return Y_eq

#calculates the variation of entropy with x, where x = m_DM/T
def dsdx(x,mDM): #mDM is the mass of the DM particle
    g_star = degrees_of_freedom(x, mDM) #g_* is the number of relativistic degrees of freedom
    ds_dx = (-6 * math.pi**2 * g_star * mDM**3)/(45 * x**4)
    return ds_dx