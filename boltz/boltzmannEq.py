#this module creates the boltzmann equation for each component

import numpy as np
import thermal.equilibriumDensities as eqDensitities
from tools.modelData import ModelData
from typing import List, Dict
from numpy.typing import ArrayLike
from tools.logger import logger

def computeDecayTerms(x: float, Y : List[float], model : ModelData) -> ArrayLike:

    mDM = model.mDM
    compDict = model.componentsDict
    mDM = compDict[model.dmPDG].mass
    T = mDM/x
    
    # The zero component is the SM, so we set dY = 0 and Y = Yeq always    
    s = eqDensitities.S(T) #entropy density at temperature T

    # Pre-compute relevant decay terms:
    Dij = np.zeros((len(Y),len(Y)))
    for comp_i in compDict.values():
        if comp_i.ID == 0: # skip SM
            continue
        if comp_i.decays is None:
            continue
        if comp_i.totalwidth <= 0.0:
            continue
        i = comp_i.ID
        Y_i = Y[i]
        Yeq_i = comp_i.Yeq(T)
        gamma_i = comp_i.gammaInv(T)
        if np.isnan(gamma_i) or gamma_i > 1.0:
            logger.error(f'Error computing gamma for {comp_i}: T = {T}, x = {comp_i.mass/T}, gamma = {gamma_i}')
        #loop over all the possible decays
        for daughter_pdgs,br in comp_i.decays.items():
            daughters = [compDict[pdg] for pdg in daughter_pdgs]             
            daughter_ids = [daughter.ID for daughter in daughters]
            for j in sorted(np.unique(daughter_ids)):
                Yprod = np.prod([Y[daughter.ID]/daughter.Yeq(T) for daughter in daughters])  
                Dij[i,j] += gamma_i*comp_i.totalwidth*br*(Y_i - Yeq_i*Yprod)/s

    return Dij

def computeCollisionTerms(x: float, Y : List[float], model : ModelData) -> List[Dict[str,float]]:

    mDM = model.mDM
    compDict = model.componentsDict
    mDM = compDict[model.dmPDG].mass
    T = mDM/x
            

    coll_terms = [{} for _ in compDict.values()]
    for comp_i in compDict.values():        
        i = comp_i.ID
        if i == 0:
            continue # Skip SM
        ### Collision term (i + a -> b + c)
        for process in model.collisionProcesses:
            if not process.changesPDG(comp_i.PDG):
                continue
            sigma = process.sigmaV(x)
            if sigma <= 0.0:
                continue
            initialPDGs = process.initialPDGs
            finalPDGs = process.finalPDGs
            # The impact on the component i is proportional
            # to the number of particles "i" destroyed in the initial state
            # minus the number of particles "i" created in the final state:
            factor = initialPDGs.count(comp_i.PDG)-finalPDGs.count(comp_i.PDG)
            # Compute collision term:
            a_pdg,b_pdg = initialPDGs
            c_pdg, d_pdg = finalPDGs
            a = compDict[a_pdg]
            b = compDict[b_pdg]
            c = compDict[c_pdg]
            d = compDict[d_pdg]    
            r_eq = a.Yeq(T)*b.Yeq(T)
            if r_eq > 0.0:
                r_eq = r_eq/(c.Yeq(T)*d.Yeq(T))
            # Multiply sigma by Yeq_i*Yeq_j for convenience
            # then we just need to multiply by ratios
            Cabcd = factor*sigma*(Y[a.ID]*Y[b.ID] - Y[c.ID]*Y[d.ID]*r_eq)
            coll_terms[i].setdefault(process.name,0.0)
            coll_terms[i][process.name] -= Cabcd

    return coll_terms

def dYdx(x: float, Y : List[float], model : ModelData):
    """
    Boltzmann equations for the BSM components. Assumes the energy density is dominated by radiation
    and the SM degrees of freedom.

    :param x: Evolution variable, x = mDM/T
    :param Y: List of yields for each BSM componnent (zero component is the SM)
    :param model: Model data object holding information about the particles and collision processes.
    """

    mDM = model.mDM
    compDict = model.componentsDict
    mDM = compDict[model.dmPDG].mass
    T = mDM/x
    # The zero component is the SM, so we set dY = 0 and Y = Yeq always
    Y[0] = compDict[0].Yeq(T)
    H = eqDensitities.H(T) #hubble rate at temperature T
    dsdx = np.abs(eqDensitities.dSdx(x, mDM)) #variation of entropy with x

    Dij = computeDecayTerms(x,Y,model)
    coll_i = computeCollisionTerms(x,Y,model)

    #Compute derivative 
    dY = np.zeros(len(Y)) #list with all the rhs for all the components, ordered as in the comp_list
    for comp_i in compDict.values():        
        i = comp_i.ID
        if i == 0:
            continue # Skip SM
        ## Decay and inverse decay terms ( i <-> j + ...)
        dec_term = -np.sum(Dij[i,:])
        ## Injection and inverse injection terms ( j <-> i + ...)
        inj_term = np.sum(Dij[:,i])
        ### Collision term (i + a <-> b + c)
        coll_term = np.sum(list(coll_i[i].values()))
        ### Divide by expansion rate
        dec_term = (dec_term/(3*H))*dsdx
        inj_term = (inj_term/(3*H))*dsdx
        coll_term = (coll_term/(3*H))*dsdx
        logger.debug(f'{comp_i} : x = {x:1.2e}, dsdx/3H = {dsdx/(3*H):1.2e}, Decay = {dec_term:1.3e}, Injection = {inj_term:1.3e} and Collision = {coll_term:1.3e}')
        dY[i] = (dec_term + inj_term + coll_term)
    dYStr = ','.join([f"{dy:1.2e}" for dy in dY])
    YStr = ','.join([f"{y:1.2e}" for y in Y])
    logger.debug(f'Result: x = {x:1.2e}, Y = {YStr}, dY = {dYStr}\n')
    return dY
