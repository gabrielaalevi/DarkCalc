from tools.modelData import ModelData
from boltz.boltzmannEq import dYdx, computeCollisionTerms,computeDecayTerms
from scipy.integrate import solve_ivp
import numpy as np
from tools.configParserWrapper import ConfigParserExt
import os
from tools.logger import setLogLevel

parfile = './input_parameters.ini'
parser = ConfigParserExt(inline_comment_prefixes="#")   
ret = parser.read(parfile)
parser = parser.expandLoops()[0]
parserDict = parser.toDict(raw=False)

pars = parser['SolverParameters']
atol = 0.0
rtol = 1e-4
T0 = pars['T0']
Tf = pars['Tf']
method = 'Radau'
nsteps = 500
outputFolder = os.path.abspath(parser['Options']['outputFolder'])
bannerFile = os.path.join(outputFolder,'darkcalc_banner.txt')
dm = parser['Model']['darkmatter']
bsmList = []
if 'bsmParticles' in parser['Model']:
    bsmList = str(parser['Model']['bsmParticles']).split(',')
model = ModelData(dmPDG=dm, bsmPDGList=bsmList, bannerFile=bannerFile)
    

compDict = model.componentsDict
mDM = compDict[model.dmPDG].mass
x0 = mDM/T0
xf = mDM/Tf  
xvals = np.geomspace(x0,xf,1000)  

# Set initila conditions
y0 = np.array([comp.Yeq(T0) for comp in compDict.values()])    
if 'initialConditions' in pars:
    # Set initila conditions
    initialCond = pars['initialConditions']    
    for label,comp_y0 in initialCond.items():
        pdg = model.convert2PDG(label)
        comp = compDict[pdg]            
        if isinstance(comp_y0,float):
            y0[comp.ID] = comp_y0
        elif comp_y0.lower() in ['eq', 'equilibrium']:
            continue # Already set
        elif comp_y0.lower() == 'zero':
            y0[comp.ID] = 1e-20

# # Run solve following Gabriela's prescription
sol = solve_ivp(dYdx, [x0,7], y0, args=(model,), atol = 10**(-12), rtol = 10**(-12), method='BDF', 
                 t_eval=np.geomspace(x0,7,50))
y = sol.y[:]
x = sol.t[:]
sol = solve_ivp(dYdx, [sol.t[-1],12], sol.y[:,-1], args=(model,), atol = 10**(-10), rtol = 10**(-10), method='BDF', max_step = 0.08,
                 t_eval=np.geomspace(7,12,50))
print(y, sol.y)
y_sol = sol.y[:]
x_sol = sol.t[:]
y = np.hstack((y,y_sol))
x = np.hstack((x,x_sol))
sol = solve_ivp(dYdx, [sol.t[-1],40], sol.y[:,-1], args=(model,), atol = 10**(-12), rtol = 10**(-12), method='BDF', max_step = 0.08,
                 t_eval=np.geomspace(12,40,50))
y_sol = sol.y[:]
x_sol = sol.t[:]
y = np.hstack((y,y_sol))
x = np.hstack((x,x_sol))
sol = solve_ivp(dYdx, [sol.t[-1],xf], sol.y[:,-1], args=(model,), atol = 10**(-10), rtol = 10**(-10), method='BDF', max_step = 0.04,
                 t_eval=np.geomspace(40,xf,50))
y_sol = sol.y[:]
x_sol = sol.t[:]
y = np.hstack((y,y_sol))
x = np.hstack((x,x_sol))

#sol = solve_ivp(dYdx, [x0,xf], y0, args=(model,), atol = atol, rtol = rtol, method=method, t_eval=xvals)
print(sol)
#x = sol.t
#y = sol.y
sol_matrix = np.zeros((len(x),len(compDict)+1))

for i in range(0, len(x)):
    sol_matrix[i][0] = x[i]
    sol_matrix[i][1] = y[0][i]
    sol_matrix[i][2] = y[1][i]
    sol_matrix[i][3] = y[2][i]

np.savetxt('CDFO result', sol_matrix)

s = 2889.2 #entropydensity(x, mDM)
Ytot = y[2][-1] #sum(y[1:,-1])
rho_c = 1.05 * 10**(-5)
omegah2 = ((s * mDM * Ytot)/rho_c)
print(omegah2)