import math
import numpy as np
import Components

#This module sets the model's main parameters, that define the DM and mediator characteristics

#parameter mDM: mass of DM particle, in GeV
#parameter gDM: degrees of freedom for the DM particle, in GeV
#parameter mMediator: mass of the mediator envolved
#parameter gMediator: degrees of freedom for the Mediator particle
#parameter decaywidthMediator: decay width of the Mediator, in GeV
#parameter MediatorDecayProducts: list holding all the possible products of the decay of the Mediator, written with the IDs
#parameter MediatorDecayBR: list holding all the Branching Ratios (as functions of T) for the possible decays of the 
#Mediator, in the same order as the MediatorDecayProducts
#parameter DMcollision: list holding all the possible collision partners for the DM particle
#parameter DMcollisionProducts: list holding all the products of the possible collisions the DM particle can participate.
#must be in the same order as the DMcollision parameter
#parameter DMcollisionSigmaV: list of functions holding all the thermally averaged cross-sections for the collisions the
#DM particle can participate. must be in the same order as the DMcollision and the DMcollisionProducts parameters
#parameter Mediatorcollision: list holding all the possible collision partners for the Mediator particle
#parameter MediatorcollisionProducts: list holding all the products of the possible collisions the Mediator particle can participate.
#must be in the same order as the Mediatorcollision parameter
#parameter MediatorcollisionSigmaV: list of functions holding all the thermally averaged cross-sections for the collisions the
#Mediator particle can participate. must be in the same order as the Mediatorcollision and the MediatorcollisionProducts parameters
#parameter T: interval of temperature for integration, in GeV. Must be in decreasing order
#parameter FO: boolean. if True, DM starts in thermal equilibrium with the plasma (a freeze-out model). If False, DM
#population is virtually zero in early times (freeze-in model).

T = np.linspace(600, 10, 9000)
mDM = 600
gDM = 1
mMediator = 750
gMediator = 1
decaywidthMediator = lambda T: 0
DecayProductsMediator = [[0,2], [0, 0], [2, 2]]
DecayBRMediator = (lambda T: 0, lambda T:0, lambda T: 0)
CollisionDM = (0, 1, 2)
CollisionProductsDM = [[1,2], [2,2], [1,2]]
CollisionSigmaVDM = (lambda T: 10e-14 * T, lambda T: 0, lambda T: 0)
CollisionMediator = (1, 0, 2)
CollisionProductsMediator = [[2, 2], [2,2], [0, 2]]
CollisionSigmaVMediator = (lambda T: 0, lambda T: 0, lambda T: 0)
ActiveList = [True, True]
FO = True