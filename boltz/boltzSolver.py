
from scipy.integrate import solve_ivp,odeint
from boltz.boltzmannEq import dYdx
from typing import List
from tools.modelData import ModelData
import numpy as np
from scipy.integrate import OdeSolution
from tools.logger import logger


def runSolver(parser : dict, model : ModelData) -> OdeSolution:
    
    logger.debug(f'Solving Boltzmann equations for {model}')

    pars = parser['SolverParameters']
    atol = pars['atol']
    rtol = pars['rtol']
    T0 = pars['T0']
    Tf = pars['Tf']
    method = pars['method']
    nsteps = pars['nsteps']
    compDict = model.componentsDict
    mDM = compDict[model.dmPDG].mass
    x0 = mDM/T0
    xf = mDM/Tf    

    
    # Initialize all components in equilibrium
    y0 = np.array([comp.Yeq(T0) for comp in compDict.values()])    
    if 'initialConditions' in pars:
        # Set initila conditions
        initialCond = pars['initialConditions']    
        for label,comp_y0 in initialCond.items():
            pdg = model.convert2PDG(label)
            comp = compDict[pdg]            
            if isinstance(comp_y0,float):
                y0[comp.ID] = comp_y0
            elif comp_y0.lower() in ['eq', 'equilibrium', 'thermal']:
                continue # Already set
            elif type(comp_y0) == int or type(comp_y0) == float:
                y0[comp.ID] = comp_y0
            else:
                logger.error(f"Could not set initial condition to {comp_y0}")
                return False
        
    xvals = np.geomspace(x0,xf,nsteps)

    solution = solveBoltzEqs(xvals,Y0=y0,model=model,
                             method=method,atol=atol,rtol=rtol)
    
    return solution

def solveBoltzEqs(xvals : List[float], Y0 : List[float], model : ModelData,
                  method : str = 'Radau', 
                  atol : float = 1e-10,
                  rtol : float = 1e-10,):


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    nsteps = len(xvals)
    #solving the Boltzmann equation

    sol = solve_ivp(dYdx, [x0,7], Y0, args=(model,), atol = 10**(-12), rtol = 10**(-12), method=method, max_step = 0.08,
                 t_eval=np.geomspace(x0,7,nsteps/4))
    y = sol.y[:]
    x = sol.t[:]
    if sol.success:
        sol2 = solve_ivp(dYdx, [sol.t[-1],12], sol.y[:,-1], args=(model,), atol = 10**(-10), rtol = 10**(-10), method=method, max_step = 0.08,
                 t_eval=np.geomspace(sol.t[-1],12,nsteps/4))
        y_sol = sol2.y[:]
        x_sol = sol2.t[:]
        y = np.hstack((y,y_sol))
        x = np.hstack((x,x_sol))
        if sol2.success:
            sol3 = solve_ivp(dYdx, [sol2.t[-1],110], sol2.y[:,-1], args=(model,), atol = 10**(-13), rtol = 10**(-13), method=method, max_step = 0.04,
                 t_eval=np.geomspace(sol2.t[-1],110,nsteps/4))
            y_sol = sol3.y[:]
            x_sol = sol3.t[:]
            y = np.hstack((y,y_sol))
            x = np.hstack((x,x_sol))
            if sol3.success:
                sol4 = solve_ivp(dYdx, [sol3.t[-1],xf], sol3.y[:,-1], args=(model,), atol = atol, rtol = rtol, method=method, max_step = 0.08,
                        t_eval=np.geomspace(sol3.t[-1],xf,nsteps/4))
                y_sol = sol4.y[:]
                x_sol = sol4.t[:]
                y = np.hstack((y,y_sol))
                x = np.hstack((x,x_sol))
                return sol4, x, y
            else:
                return sol3, x, y
        else:
            return sol2, x, y
    else:
        return sol, x, y