#solving the ODE

from scipy.integrate import solve_ivp
import numpy as np
from boltzmannEq import boltz
from modelParameters import comp_names, x, mDM, SM
import numpy as np

#creating the y0 array, that holds the initial value of yield for each particle, and loop to fill it
y0 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    if comp.in_equilibrium == 1:
        y0[i] = comp.equilibriumyield(x[0], mDM)
    if comp.in_equilibrium == 0:
        y0[i] = 0.0001
#solving the Boltzmann equation
sol1 = solve_ivp(boltz, [x[0], 9.0], y0, args=(comp_names, SM, mDM), atol = 10**(-12), rtol = 10**(-12), method='BDF')
y02 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y02[i] = sol1.y[i][-1]
sol2 = solve_ivp(boltz, [9.0, 10.5],y02, args=(comp_names, SM, mDM), atol = 10**(-12), rtol = 10**(-12), method='LSODA', max_step = 0.08)
y03 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y03[i] = sol2.y[i][-1]
sol3 = solve_ivp(boltz, [10.5, 28], y03, args=(comp_names, SM, mDM), atol = 10**(-12), rtol = 10**(-12), method='LSODA', max_step = 0.08)
y04 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y04[i] = sol3.y[i][-1]
sol4 = solve_ivp(boltz, [28, x[-1]],y04, args=(comp_names, SM, mDM), atol = 10**(-12), rtol = 10**(-12), method='LSODA', max_step = 0.08)
solt = [*sol1.t, *sol2.t, *sol3.t, *sol4.t]
#in here, we create the header for the final txt file
header = ' '
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    header = header + '  Y' + comp.label + '  Y_eq' + comp.label
#calculating all the equilibrium yields for the participating particles
y_eq= np.zeros((len(comp_names), len(solt)))
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    for j in range(0, len(solt)):
        y_eq[i][j] = comp.equilibriumyield(solt[j], mDM)

