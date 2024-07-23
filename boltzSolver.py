#Solver for the Boltzmann Equation
from scipy.integrate import odeint
from boltzmannEq import boltz
from components import Component
import auxFunc
from modelParameters import particlesDict, Tlist
import numpy as np
comp_list = []
for i in range(0, len(particlesDict)):
    comp = Component(particlesDict[i]['label'], particlesDict[i]['active'], particlesDict[i]['ID'], particlesDict[i]['g'], particlesDict[i]['mass'], particlesDict[i]['decaywidth'], particlesDict[i]['decayreactions'], particlesDict[i]['collisions'], particlesDict[i]['equilibrium'])
    comp_list.append(comp)
y0 = [comp.equilibriumyield(Tlist[0]) for comp in comp_list]
for k in range(0, len(comp_list)):
    comp = comp_list[k]
    if comp.label == 'DM':
        x = comp.mass/Tlist
sol = odeint(boltz, y0, x, args=(comp_list,), atol = 10**(-13), rtol = 10**(-14))
