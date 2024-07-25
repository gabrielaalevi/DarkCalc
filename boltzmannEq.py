#this module creates the boltzmann equation for each component

import numpy as np
from scipy.special import kn
import auxFunc

def boltz(Y, x, comp_list):
    dY = np.zeros(len(comp_list)) #list with all the rhs for all the components, ordered as in the comp_list
    for i in range(0, len(comp_list)):
        comp = comp_list[i]
        if comp.label == 'DM':
                mDM = comp.mass
    H = auxFunc.hubblerate(x, mDM) #hubble rate at temperature T
    s = auxFunc.entropydensity(x, mDM) #entropy density at temperature T
    dsdx = auxFunc.dsdx(x, mDM) #variation of entropy with x
    T = mDM/x
    for i in range(0, len(comp_list)):
    #loop parsing through all components
        comp = comp_list[i]
        #decay term for component i
        if comp.decayreactions != 0:
            for j in range(0, len(comp.decayreactions)):
            #loop parsing through all the possible decays of component i
                dec_term = 1
                product = []
                for l in range(0, len(comp.decayreactions[j]['products'])):
                #loop parsing through the products of decay j
                    for k in range(0, len(comp_list)):
                        #loop parsing through all the components, to see which ones are products in the decay
                        comp_product = comp_list[k]
                        if (comp.decayreactions[j]['products'][l]) == comp_product.label:
                            dec_term *= Y[comp_product.ID]/comp_product.equilibriumyield(x, mDM)
                            product.append(comp_product.ID)
                dY[i] += - (1/(3 * H)) * dsdx * (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM) * comp.decayreactions[j]['br'] *dec_term)
                for l in range(0, len(product)):
                    #adding the correspondent source term to the dY equation of each product
                    dY[(product[l])] += (1/(3 * H)) * dsdx * (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM) * comp.decayreactions[j]['br'] * dec_term)
        #collision term for component i
        for m in range(0, len(comp.collisions)):
        #loop parsing through all the collisions component i can participate in as a reagent
            col_term = 1
            product = []
            partner = comp.collisions[m]['partner']
            for n in range(0, len(comp.collisions[m]['products'])):
            #loop parsing through the products of the collision m
                for p in range(0, len(comp_list)):
                #loop parsing through all the components, to see which ones are the products of the collision
                    comp_product = comp_list[p]
                    if partner == comp_product.label:
                    #if the component we are analysing is the partner, we save it as the partner
                        partner = comp_product
                    if (comp.collisions[m]['products'][n]) == comp_product.label:
                    #if the component we are analysing is one of the products of the reaction, we add to the col_term and save its ID
                        col_term *= Y[comp_product.ID]/comp_product.equilibriumyield(x, mDM)
                        product.append(comp_product.ID)
                if partner.label != 'SM':
                    dY[i] += - (1/(3 * H)) * dsdx* comp.collisions[m]['sigmaV'](T) *((Y[i] * Y[partner.ID])/(comp.equilibriumyield (x, mDM) * partner.equilibriumyield(x, mDM)) - col_term)
                    dY[partner.ID] += - (1/(3*H)) * dsdx * comp.collisions[m]['sigmaV'](T) * ((Y[i] * Y[partner.ID])/(comp.equilibriumyield(x, mDM) * partner.equilibriumyield(x, mDM)) - col_term)
                    for p in range(0, len(product)):
                        dY[(product[p])] += (1/(3 * H)) * dsdx * comp.collisions[m]['sigmaV'](T) * ((Y[i]*Y[partner.ID])/(comp.equilibriumyield(x, mDM) * partner.equilibriumyield(x, mDM)) - col_term)
                if partner.label == 'SM':
                    dY[i] += (1/(3*H)) * dsdx * comp.collisions[m]['sigmaV'](T) * ((Y[i]/ comp.equilibriumyield(x, mDM)) - col_term)
                    for p in range(0, len(product)):
                        dY[(product[p])] += (1/(3 * H)) * dsdx * comp.collisions[m]['sigmaV'](T) * ((Y[i]/comp.equilibriumyield(x, mDM)) - col_term)
    return dY
                
            
        
