#example main

import math
import numpy as np
import BoltzmannEq
import Components
import AuxFunc
import ModelParameters
from scipy.integrate import odeint
import matplotlib.pyplot as plt


dm = Components.Component(label = 'dm', ID = 0, active = ModelParameters.ActiveList[0], g = ModelParameters.gDM, mass = ModelParameters.mDM, decaywidth = 0, decayproducts = 0, br = 0, collision = ModelParameters.CollisionDM, collisionproducts = ModelParameters.CollisionProductsDM, sigmaV = ModelParameters.CollisionSigmaVDM)
mediator = Components.Component(label = 'mediator', ID = 1, active = ModelParameters.ActiveList[1], g = ModelParameters.gMediator, mass = ModelParameters.mMediator, decaywidth = ModelParameters.decaywidthMediator, decayproducts = ModelParameters.DecayProductsMediator, br = ModelParameters.DecayBRMediator, collision = ModelParameters.CollisionMediator, collisionproducts = ModelParameters.CollisionProductsMediator, sigmaV = ModelParameters.CollisionSigmaVMediator)

comp_list = [dm, mediator]

y0_dm = AuxFunc.equilibrium_yield(ModelParameters.mDM, ModelParameters.T[0], ModelParameters.gDM)
y0_mediator = AuxFunc.equilibrium_yield(ModelParameters.mMediator, ModelParameters.T[0], ModelParameters.gMediator)
x = ModelParameters.mDM/ModelParameters.T

sol = odeint(BoltzmannEq.BoltzmannEq.Boltz, [y0_dm,y0_mediator], x, args=(comp_list,), atol = 10**(-16), rtol = 10**(-14))

sol_dm = []
sol_mediator = []

for i in range(0, len(ModelParameters.T)):
    sol_dm.append(sol[i][0])
    sol_mediator.append(sol[i][1])

plt.plot(x, sol_dm)
plt.yscale('log')
plt.xscale('log')
plt.grid()