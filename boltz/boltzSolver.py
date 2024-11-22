
from scipy.integrate import solve_ivp, odeint
from boltz.boltzmannEq import boltz
from typing import List
from modelData import ModelData


def solveBoltzEqs(xvals : List[float], Y0 : List[float], model : ModelData,
                  method : str = 'BDF', 
                  atol : float = 1e-12,
                  rtol : float = 1e-12,):


    # Initial conditions
    x0, xf = xvals[0],xvals[-1]
    #solving the Boltzmann equation
    sol = solve_ivp(boltz, [x0,xf], Y0, args=(model,), atol = atol, rtol = rtol, method=method)
    # sol1= odeint(boltz, y0, xvals, args=(model,), tfirst=True)

    return sol