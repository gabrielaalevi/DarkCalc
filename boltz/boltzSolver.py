
from scipy.integrate import solve_ivp,odeint
from boltz.boltzmannEq import dYdx
from scipy.integrate import OdeSolution
from typing import List,Tuple,Union
from tools.modelData import ModelData
import numpy as np
from numpy.typing import NDArray
from tools.logger import logger
import warnings


def runSolver(parser : dict, model : ModelData) -> Tuple[Union[OdeSolution,None],NDArray,NDArray]:
    
    logger.debug(f'Solving Boltzmann equations for {model}')

    pars = parser['SolverParameters']
    atol = pars.get('atol',1e-10)
    rtol = pars.get('rtol',1e-4)
    T0 = pars['T0']
    Tf = pars['Tf']
    method = pars.get('method','Radau')
    nsteps = pars.get('nsteps',10000)
    max_step = pars.get('max_step',0.1)
    compDict = model.componentsDict
    mDM = compDict[model.dmPDG].mass
    x0 = mDM/T0
    xf = mDM/Tf

    mass_max = max([ptc.mass for ptc in compDict.values()])
    with warnings.catch_warnings():
        warnings.filterwarnings('error', category=RuntimeWarning) # Treat RuntimeWarning as an error
        try:
            np.exp(mass_max/Tf)
        except RuntimeWarning:
            msg = f"Final temperature value {Tf:1.2e} is well below"
            msg += f"the largest BSM mass ({mass_max:1.3f}) and may cause numerical instabilities"
            logger.warning(msg)
    
    
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
                             method=method,atol=atol,rtol=rtol,
                             max_step=max_step)
    
    return solution

def solveBoltzEqs(xvals : NDArray, Y0 : NDArray, 
                  model : ModelData,
                  method : str = 'Radau', 
                  atol : float = 1e-10,
                  rtol : float = 1e-10,
                  max_step : float = 0.1,) -> Tuple[Union[OdeSolution,None],NDArray,NDArray]:


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    nsteps = len(xvals)
    n_eval = int(nsteps/4)
    #solving the Boltzmann equation

    # It is usually convenient to break the solution into smaller intervals:
    xfList = [7.,12.,110.,xf]
    xfList = [x for x in xfList[:] if x0 < x <= xf]


    x = np.array([])
    y = np.array([[] for _ in Y0])
    sol = None
    y0=Y0[:]
    for xf_val in xfList:
        logger.debug(f'Solving in the interval {x0} to {xf_val}')
        sol = solve_ivp(dYdx, [x0,xf_val], y0=y0, args=(model,), 
                        atol = atol, rtol = rtol, first_step=max_step/5,
                        method=method, max_step = max_step,
                        t_eval=np.geomspace(x0,xf_val,n_eval))
        if not sol.success:
            return sol, x, y
        
        y = np.hstack((y,sol.y[:]))
        x = np.hstack((x,sol.t[:]))
        x0 = sol.t[-1]
        y0 = sol.y[:,-1]
    
    return sol, x, y

def solveBoltzEqs_old(xvals : List[float], Y0 : List[float], model : ModelData,
                  method : str = 'Radau', 
                  atol : float = 1e-10,
                  rtol : float = 1e-10,
                  max_step : float = 0.1,):


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    nsteps = len(xvals)
    n_eval = int(nsteps/4)
    #solving the Boltzmann equation

    sol = solve_ivp(dYdx, [x0,7], Y0, args=(model,), atol = 10**(-12), rtol = 10**(-12), method=method, max_step = 0.08,
                 t_eval=np.geomspace(x0,7,n_eval))
    y = sol.y[:]
    x = sol.t[:]
    if sol.success:
        sol2 = solve_ivp(dYdx, [sol.t[-1],12], sol.y[:,-1], args=(model,), atol = 10**(-10), rtol = 10**(-10), method=method, max_step = 0.08,
                 t_eval=np.geomspace(sol.t[-1],12,n_eval))
        y_sol = sol2.y[:]
        x_sol = sol2.t[:]
        y = np.hstack((y,y_sol))
        x = np.hstack((x,x_sol))
        if sol2.success:
            sol3 = solve_ivp(dYdx, [sol2.t[-1],110], sol2.y[:,-1], args=(model,), atol = 10**(-13), rtol = 10**(-13), method=method, max_step = 0.04,
                 t_eval=np.geomspace(sol2.t[-1],110,n_eval))
            y_sol = sol3.y[:]
            x_sol = sol3.t[:]
            y = np.hstack((y,y_sol))
            x = np.hstack((x,x_sol))
            if sol3.success:
                sol4 = solve_ivp(dYdx, [sol3.t[-1],xf], sol3.y[:,-1], args=(model,), atol = atol, rtol = rtol, method=method, max_step = 0.08,
                        t_eval=np.geomspace(sol3.t[-1],xf,n_eval))
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
    