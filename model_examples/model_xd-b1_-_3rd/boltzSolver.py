#solving the ODE

from scipy.integrate import solve_ivp
import numpy as np
from boltzmannEq import boltz
from modelParameters import comp_names, x, mDM, SM
import numpy as np
import math

#creating the y0 array, that holds the initial value of yield for each particle, and loop to fill it
y0 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    if comp.in_equilibrium == 1:
        y0[i] = comp.equilibriumyield(x[0], mDM)
    if comp.in_equilibrium == 0:
        y0[i] = 0.0001
sol1 = solve_ivp(boltz, [x[0], 9.0], y0, args=(comp_names, SM, mDM), atol = 10**(-13), rtol = 10**(-14), method='BDF')
tolv = [10**(-10), 10**(-13)] #atol and rtol
sol2 = solve_ivp(boltz, [9.0, 10.5], [sol1.y[0][-1], sol1.y[1][-1]], args=(comp_names, SM, mDM), atol = tolv[0], rtol = tolv[1], method='LSODA', max_step = 10**(-2))
sol3 = solve_ivp(boltz, [10.5, 28], [sol2.y[0][-1], sol2.y[1][-1]], args=(comp_names, SM, mDM), atol = 10**(-11), rtol = 10**(-13), method='LSODA', max_step = 10**(-2))
sol4 = solve_ivp(boltz, [28, x[-1]], [sol3.y[0][-1], sol3.y[1][-1]], args=(comp_names, SM, mDM), atol = 10**(-11), rtol = 10**(-10), max_step = 10**(-2))
# if sol2.status != 0:
#     for i in range (0,10):
#         if sol2.status !=0:
#             for k in range(1, 10):
#                 if sol2.status != 0 :
#                     rtolv = 10**(-14) * k * 10**(i)
#                     sol2 = solve_ivp(boltz, [9.3, x[-1]], [sol1.y[0][-1], sol1.y[1][-1]], args=(comp_names, SM), rtol = rtolv, atol = 10**(-14), method='LSODA')
#                     if sol2.status != 0:
#                         for j in range(0, 4):
#                             if sol2.status !=0:
#                                 for l in range(1, 10):
#                                     if sol2.status != 0:
#                                         atolv = 10**(-14) * l * 10**(j)
#                                         print('atol', atolv, 'rtol', rtolv)
#                                         sol2 = solve_ivp(boltz, [9.3, x[-1]], [sol1.y[0][-1], sol1.y[1][-1]], args=(comp_names, SM), atol = atolv, rtol = rtolv, method='LSODA')
header = '  '
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    header = header + '  Y' + comp.label + '  Y_eq' + comp.label
print(header)
