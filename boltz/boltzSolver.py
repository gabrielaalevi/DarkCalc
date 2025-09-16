
from scipy.integrate import solve_ivp,odeint
from boltz.boltzmannEq import dYdx
from scipy.integrate import OdeSolution
from typing import List,Tuple,Union
from tools.modelData import ModelData
import numpy as np
from numpy.typing import NDArray
from tools.logger import logger


def runSolver(parser : dict, model : ModelData) -> Tuple[Union[OdeSolution,None],NDArray,NDArray]:
    
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
                raise ValueError()
        
    xvals = np.geomspace(x0,xf,nsteps)

    solution = solveBoltzEqs(xvals,Y0=y0,model=model,
                             method=method,atol=atol,rtol=rtol)
    
    return solution

def solveBoltzEqs(xvals : NDArray, Y0 : NDArray, 
                  model : ModelData,
                  method : str = 'Radau', 
                  atol : float = 1e-10,
                  rtol : float = 1e-10,) -> Tuple[Union[OdeSolution,None],NDArray,NDArray]:


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    nsteps = len(xvals)
    n_eval = int(nsteps/4)
    #solving the Boltzmann equation

    # It is usually convenient to break the solution into smaller intervals:
    xfList = [7.,12.,110.,xf]
    xfList = [x for x in xfList[:] if x <= xf]


    x = np.array([])
    y = np.array([[] for _ in Y0])
    sol = None
    y0=Y0[:]
    for xf_val in xfList:
        sol = solve_ivp(dYdx, [x0,xf_val], y0=y0, args=(model,), 
                        atol = atol, rtol = rtol, 
                        method=method, max_step = 0.1,
                        t_eval=np.geomspace(x0,xf_val,n_eval))
        if not sol.success:
            return sol, x, y
        
        y = np.hstack((y,sol.y[:]))
        x = np.hstack((x,sol.t[:]))
        x0 = sol.t[-1]
        y0 = sol.y[:,-1]
    
    return sol, x, y
    