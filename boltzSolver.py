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
sol = solve_ivp(boltz, [x[0], x[-1]], y0, args=(comp_names, SM), atol=10**(-11), rtol=10**(-12), method='LSODA')