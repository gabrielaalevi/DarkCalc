#Solver for the Boltzmann Equation
from scipy.integrate import solve_ivp
import numpy as np
from boltzmannEq import boltz
from components import Component
import auxFunc
from modelParameters import comp_names, Tvalues, x
from inputParameters import SM
import numpy as np

#loop to find the particle with the 'DM' label, and calculate x
for k in range(0, len(comp_names)):
    comp = comp_names[k]
    if comp.label == 'DM':
        x = comp.mass/Tvalues
        mDM = comp.mass

#creating the y0 array, that holds the initial value of yield for each particle, and loop to fill it
y0 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    if comp.in_equilibrium == 1:
        y0[i] = comp.equilibriumyield(x[0], mDM)
    if comp.in_equilibrium == 0:
        y0[i] = 0.0001

#solving the ODE
sol1 = solve_ivp(boltz, [x[0], 9.5],  y0, args=(comp_names, SM), method='BDF', dense_output= True, rtol = 10**(-14), atol = 10**(-14))
sol2 = solve_ivp(boltz, [9.5, 10.5], (sol1.y[0][-1], sol1.y[1][-1]), args=(comp_names, SM), method='BDF', dense_output= True, rtol = 10**(-14), atol = 10**(-14))
sol3 = solve_ivp(boltz, [10.5, 15], (sol2.y[0][-1], sol2.y[1][-1]), args=(comp_names, SM), method='BDF', dense_output= True, rtol = 10**(-14), atol = 10**(-14))
sol4 = solve_ivp(boltz, [15, x[-1]], (sol3.y[0][-1], sol3.y[1][-1]), args=(comp_names, SM), method='BDF', dense_output= True, rtol = 10**(-12), atol = 10**(-14))
#np.savetxt('sol.csv', sol)
