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


def degrees_of_freedom(T): #temperature in GeV
    if T >= 0.3:
        g = 100
    elif T >= 0.1:
        g = 10
    elif 0.1 > (T):
        g = 3
    return g


# In[4]:


def hubble(T): #T for temperature, g for the relativistic number of degrees of freedom
    g = degrees_of_freedom(T)
    hubble = ((np.pi**2)/90)**(1/2) * g**(1/2) * T**2 / m_planck #hubble constant, varying with temperature
    return hubble


# In[5]:


def entropy(T): #entropy density, varying with temperature T and degrees of freedom g
    g = degrees_of_freedom(T)
    s = (2/45) * (np.pi**2) * g * T**3
    return s


# In[6]:


def equilibrium_yield(m, T, g,b):
    #m represents particle's mass, T is temperature in GeV, g represents, particles degrees of freedom, 
    #and s represents entropy density. b is a parameter, b=0 indicates a fermion, and b=1 indicates a bóson
    if T >= 0.3:
        g_SM = 100
    elif T >= 0.1:
        g_SM = 10
    elif 0.1 > (T):
        g_SM = 3
    s = (2/45) * (np.pi**2) * g_SM * T**3
    if (m/T)>10:
        Y_eq = (g * (m * T/ (2 * np.pi))**(3/2) * np.exp(-m/T))/s
    elif (m/T)>(2/3):
        Y_eq = (45 * g * (m/T)**2 * kn(2,(m/T)))/(4 * np.pi**4 * g_SM)
    elif (m/T)<=(2/3):
        if b == 0:
            Y_eq = ((3/4) * (zeta(3)/np.pi**2) * g * (T)**3)/s
        if b==1:
            Y_eq = ((zeta(3)/np.pi**2) * g * (T)**3)/s
    return Y_eq


# In[7]:


def dsdx(x,m):
    g = degrees_of_freedom(m/x)
    ds_dx = (-6 * np.pi**2 * g * m**3)/(45 * x**4)
    return ds_dx

