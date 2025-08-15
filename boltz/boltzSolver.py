
from scipy.integrate import solve_ivp
from boltz.boltzmannEq import dYdx
from tools.modelData import ModelData
import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import OdeSolution
from tools.logger import logger
from typing import Union,List


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
            elif comp_y0.lower() in ['eq', 'equilibrium']:
                continue # Already set
            elif comp_y0.lower() == 'zero':
                y0[comp.ID] = 1e-20
            else:
                raise ValueError(f"Could not set initial condition to {comp_y0}")
        
    xvals = np.geomspace(x0,xf,nsteps)

    solution = solveBoltzEqs(xvals,Y0=y0,model=model,
                             method=method,atol=atol,rtol=rtol)
    
    return solution

def solveBoltzEqs(xvals : List, Y0 : Union[List,ArrayLike], 
                  model : ModelData,
                  method : str = 'Radau', 
                  atol : float = 0.0,
                  rtol : float = 1e-3,):


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    #solving the Boltzmann equation
    sol = solve_ivp(dYdx, [x0,xf], Y0, args=(model,), atol = atol, 
                    rtol = rtol, method=method, t_eval=xvals)

    return sol