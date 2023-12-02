#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
from scipy.special import zeta, kn


# In[2]:


#constants
m_planck = 2.4 * 10**(18) #reduced planck mass, in GeV


# In[3]:


#calculates the number of relativistic degrees of freedom for temperature T (in GeV), using solely the SM degrees of freedom
def degrees_of_freedom(T):
    if T >= 0.3:
        g_* = 100
    elif T >= 0.1:
        g_* = 10
    elif 0.1 > (T):
        g_* = 3
    return g_*


# In[4]:


#calculates the hubble rate H based on the temperature T (in GeV)
def hubble(T):
    g_* = degrees_of_freedom(T) #g_* is the number of relativistic degrees of freedom
    hubble = ((np.pi**2)/90)**(1/2) * (g_*)**(1/2) * T**2 / m_planck #hubble constant, varying with temperature
    return hubble


# In[5]:


#calculates the entropy density s, depending on the temperature T (in GeV)
def entropy(T):
    g_* = degrees_of_freedom(T) #g_* is the number of relativistic degrees of freedom
    s = (2/45) * (np.pi**2) * (g_*) * T**3 #entropy density
    return s


# In[6]:


#calculates the equilibrium yield for a certain particle, using Fermi-Dirac or Bose-Einstein Statistics.
def equilibrium_yield(m, T, g,b):
    #m represents particle's mass, T is temperature in GeV, g represents the particle's degrees of freedom, 
    #and s represents entropy density. b is an integer parameter, b=0 indicates a fermion, and b=1 indicates a boson
    g_* = degrees_of_freedom(T)#g_* is the number of relativistic degrees of freedom
    s = (2/45) * (np.pi**2) * (g_*) * T**3 #entropy density
    if (m/T)>10: #non-relativistic regime
        Y_eq = (g * (m * T/ (2 * np.pi))**(3/2) * np.exp(-m/T))/s
    elif (m/T)>(2/3): #semi-relativistic regime
        Y_eq = (45 * g * (m/T)**2 * kn(2,(m/T)))/(4 * np.pi**4 * g_*)
    elif (m/T)<=(2/3): #relativistic regime
        if b == 0:
            Y_eq = ((3/4) * (zeta(3)/np.pi**2) * g * (T)**3)/s
        if b==1:
            Y_eq = ((zeta(3)/np.pi**2) * g * (T)**3)/s
    return Y_eq


# In[7]:


#calculates the variation of entropy with x, where x = m_DM/T
def dsdx(x,m): #m is the mass of the DM particle
    g_* = degrees_of_freedom(m/x) #g_* is the number of relativistic degrees of freedom
    ds_dx = (-6 * np.pi**2 * g_* * m**3)/(45 * x**4)
    return ds_dx

