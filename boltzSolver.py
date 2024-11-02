#solving the ODE

from scipy.integrate import solve_ivp
import numpy as np
from boltzmannEq import boltz, debug_func
from modelParameters import comp_names, x, mDM, SM, debug_version, name_file
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
sol1= solve_ivp(boltz, [x[0], 9], y0, args=(comp_names, SM, mDM, x), atol = 10**(-12), rtol = 10**(-12), method='BDF')
y02 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y02[i] = sol1.y[i][-1]
sol2= solve_ivp(boltz, [9, 10.5],y02, args=(comp_names, SM, mDM, x), atol = 10**(-12), rtol = 10**(-12), method='LSODA', max_step = 0.08)
y03 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y03[i] = sol2.y[i][-1]
sol3 = solve_ivp(boltz, [10.5, 28], y03, args=(comp_names, SM, mDM, x), atol = 10**(-14), rtol = 3 * 10**(-14), method='LSODA', max_step = 0.05)
y04 = np.zeros(len(comp_names))
for i in range(0, len(comp_names)):
    y04[i] = sol3.y[i][-1]
sol4 = solve_ivp(boltz, [28, x[-1]],y04, args=(comp_names, SM, mDM, x), atol = 10**(-12), rtol = 10**(-12), method='LSODA', max_step = 0.08)
solt = [*sol1.t, *sol2.t, *sol3.t, *sol4.t]
Y = np.zeros(((2*len(comp_names))+1, len(solt))) #saving the output in one matrix
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    Y[0][:] = solt
    Y[2*i+1][:] = [*sol1.y[i], *sol2.y[i], *sol3.y[i], *sol4.y[i]]
    for j in range(0, len(solt)):
        Y[(2*i) + 2][j] = comp.equilibriumyield(solt[j], mDM)
#in here, we create the header for the final txt file
header = 'x'
for i in range(0, len(comp_names)):
    comp = comp_names[i]
    header = header + '         Y_' + comp.label + '           Y_eq_' + comp.label
name = name_file + ' data.csv'
Y_trans = np.transpose(Y)
np.savetxt(name, Y_trans, header=header, delimiter='   ')

if debug_version == True:
    name = name_file + ' debug.csv'
    reactions_values = list()
    header = '             x                   '
    for i in range(0, len(solt)):
        reac = debug_func(solt[i], Y_trans[i][1:], comp_names, SM, mDM, x)
        reactions_values.append(reac)
    for i in range(0, len(comp_names)):
        comp = comp_names[i]
        if comp.decayreactions == 0:
            for m in range(0, len(comp.collisions)):
                header = header + ' ' + comp.type + str(comp.collisions[m][0]) + ' -> ' 
                if not comp.collisions[m][1] or len(comp.collisions[m][1]) == 1:
                    header = header + ' SM '
                else:
                    for p in range(0, len(comp.collisions[m][1])):
                        prod = comp.collisions[m][1][p]
                        header = header + str(prod) + '  '
        else:
            for j in range(0, len(comp.decayreactions)):
                header = header + ' ' + comp.type + ' -> ' + comp.decayreactions[j][0]
                for p in range(0, len(comp.decayreactions[j][0])):
                    prod = comp.decayreactions[j][0][p]
                    header = header + str(prod)
            for m in range(0, len(comp.collisions)):
                header = header + ' ' + comp.type + comp.collisions[m][0] + ' -> ' + comp.collisions[m][1]
                if not comp.collisions[m][1] or len(comp.collisions[m][1]) == 1:
                    header = header + ' SM                       '
                else:
                    for p in range(0, len(comp.collisions[m][1])):
                        prod = comp.collisions[m][1][p]
                        header = header + str(prod)
    np.savetxt(name, np.column_stack([reactions_values]), header = header, delimiter='   ')