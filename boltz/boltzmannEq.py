#this module creates the boltzmann equation for each component

import numpy as np
from scipy.special import kn
import thermal.equilibriumDensities as eqDensitities
from components import Component, CollisionProcesses
from typing import List,Dict

def boltz(x: float, Y : List[float], compDict : Dict[int,Component], 
          mDM : float, collisions : CollisionProcesses):
    """
    Boltzmann equations for the BSM components. Assumes the energy density is dominated by radiation
    and the SM degrees of freedom.

    :param x: Evolution variable, x = mDM/T
    :param Y: List of yields for each BSM componnent
    :param compDict: Dictionary with PDG as keys and Component objects as values
    :param mDM: Dark Matter mass
    :param collisions: CollisionProcesses object holding all the relevant thermally averaged collision cross-sections
    """

    T = mDM/x
    dY = np.zeros(len(compDict)) #list with all the rhs for all the components, ordered as in the comp_list
    H = eqDensitities.H(T) #hubble rate at temperature T
    s = eqDensitities.S(T) #entropy density at temperature T
    dsdx = eqDensitities.dSdx(x, mDM) #variation of entropy with x


    # Filter the BSM components appearing in the solutions
    compDictBSM = {pdg : comp for pdg,comp in compDict.items() 
                    if (comp.ID >= 0 and comp.ID < len(Y))}
    # Pre-coompute relevant decay terms:
    Dij = np.zeros((len(compDictBSM),len(compDictBSM)))
    for comp_i in enumerate(compDictBSM.values()):
        if comp_i.decays is None:
            continue
        if comp_i.totalwidth <= 0.0:
            continue
        i = comp_i.ID
        Y_i = Y[i]
        Yeq_i = comp_i.Yeq(T)
        gamma_i = (kn(1,comp_i.mass/T)/kn(2,comp_i.mass/T))
        f_i = gamma_i*comp_i.totalwidth/s
        #loop over all the possible decays
        for daughter_pdgs,br in comp_i.decays.items():
            daughters = [compDict[pdg] for pdg in daughter_pdgs]                
            daughter_ids = [daughter.ID for daughter in daughters]
            for j in sorted(np.unique(daughter_ids)):
                Yprod = np.prod([Y[daughter.ID]/daughter.Yeq(T) for daughter in daughters])  
                Dij[i,j] += f_i*br*(Y_i - Yeq_i*Yprod)
            

        
    #loop over all components
    for comp_i in compDictBSM.values():
        i = comp_i.ID
        Y_i = Y[i]
        ## Decay and inverse decay terms ( i -> j + ...)
        dec_term = -np.sum(Dij[i,:])
        ## Injection and inverse injection terms ( j -> i + ...)
        inj_term = np.sum(Dij[:,i])
        ### STOP HERE !!!!!!!!!!!!!!!!!
        coll_term = 0.0

        dY[i] += dec_term + inj_term + coll_term
    return dY

def debug_func(x, Y, comp_names, SM, mDM, x_new):
    nofreactions = 0 #counting all the reactions in the model
    for j in range(0, len(comp_names)):
        comp = comp_names[j]
        if comp.decayreactions == 0:
            nofreactions = nofreactions + len(comp.collisions)
        else:
            nofreactions = nofreactions + len(comp.decayreactions) + len(comp.collisions)
    reaction_terms = np.zeros((nofreactions+1))
    reaction_terms[0] = x
    counter = 1 #counting all the reactions that have been saved
    H = equilibriumDensities.hubblerate(x, mDM) #hubble rate at temperature T
    s = equilibriumDensities.entropydensity(x, mDM) #entropy density at temperature T
    dsdx = equilibriumDensities.dsdx(x, mDM) #variation of entropy with x
    T = mDM/x
    tol = abs(x - x_new[0])
    index = 0
    for p in range(0, len(x_new)):
    #loop to find the value of x closest to the one being used by odeint, necessary to call the right element of the cross section and decay width arrays
        diff = abs(x - x_new[p])
        if diff < tol:
            index = p
            tol = diff
    for i in range(0, len(comp_names)):
    #loop parsing through all components
        comp = comp_names[i]
        #decay term for component i
        if comp.decayreactions != 0:
            for j in range(0, len(comp.decayreactions)):
            #loop parsing through all the possible decays of component i
                dec_term = 1
                product = []
                for l in range(0, len(comp.decayreactions[j][0])):
                #loop parsing through the products of decay j
                    for k in range(0, len(comp_names)):
                        #loop parsing through all the components, to see which ones are products in the decay
                        comp_product = comp_names[k]
                        if (comp.decayreactions[j][0][l]) == comp_product.PDG:
                            dec_term *= Y[2* comp_product.ID]/Y[2 * comp_product.ID + 1]
                            product.append(comp_product.ID)
                reaction_terms[counter] = ((1/(3 * H))*(kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[2*i] - Y[(2*i+1)] *dec_term))/H
                counter += 1
    #collision term for component i
        for m in range(0, len(comp.collisions)):
        #loop parsing through all the collisions component i can participate in as a reagent
            col_term = 1
            product = []
            partner = comp.collisions[m][0] #saving the TYPE of the partner
            sigma_v = comp.collisions[m][2][index]
            for n in range(0, len(comp_names)):
            #loop parsing through all the particles in comp_names
                comp_product = comp_names[n]
                if partner == comp_product.type:
                #if the component we are analysing is the partner, we save it as the CLASS COMPONENT
                    partner_comp = comp_product
                for p in range(0, len(comp.collisions[m][1])):
                #loop parsing through all the products of collision m
                    if (comp.collisions[m][1][p]) == comp_product.type:
                    #if the component we are analysing is one of the products of the reaction, we add to the col_term and save its ID
                        col_term *= Y[2 * comp_product.ID]/Y[2* comp_product.ID + 1]
                        product.append(comp_product.ID)
            if sigma_v != 0:
                if partner != 'SM':
                #for self-annihilation, co-annihilation, and double conversion interactions
                    reaction_terms[counter] = ((1/(3 * H)) * dsdx * sigma_v * ((Y[2* i] * Y[2 * partner_comp.ID]) - (Y[(2*i)+1] * Y[2 * partner_comp.ID + 1] * col_term)))/H
                else:
                #for conversion interaction
                    reaction_terms[counter] = ((1/(3 * H)) * dsdx * sigma_v * ((Y[2*i]) - Y[(2*i)+1] * col_term))/H
            counter += 1
    return reaction_terms
                
            
        
