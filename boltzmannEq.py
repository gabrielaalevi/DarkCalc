#this module creates the boltzmann equation for each component

import numpy as np
from scipy.special import kn
import auxFunc
from modelParameters import mDM

def boltz(x, Y, comp_names, SM):
    dY = np.zeros(len(comp_names)) #list with all the rhs for all the components, ordered as in the comp_list
    H = auxFunc.hubblerate(x, mDM) #hubble rate at temperature T
    s = auxFunc.entropydensity(x, mDM) #entropy density at temperature T
    dsdx = auxFunc.dsdx(x, mDM) #variation of entropy with x
    T = mDM/x
    for i in range(0, len(comp_names)):
    #loop parsing through all components
        comp = comp_names[i]
        x_new = comp.xarray
        tol = abs(x - x_new[0])
        index = 0
        for p in range(0, len(x_new)):
        #loop to find the value of x closest to the one being used by odeint, necessary to call the right element of the cross section and decay width arrays
            diff = abs(x - x_new[p])
            if diff < tol:
                index = p
                tol = diff
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
                            dec_term *= Y[comp_product.ID]/comp_product.equilibriumyield(x, mDM)
                            product.append(comp_product.ID)
                dY[i] += - (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM) *dec_term)
                for l in range(0, len(product)):
                #adding the correspondent source term to the dY equation of each product
                    dY[(product[l])] += (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM)* dec_term)
        
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
                    partner = comp_product
                for p in range(0, len(comp.collisions[m][1])):
                #loop parsing through all the products of collision m
                    if (comp.collisions[m][1][p]) == comp_product.type:
                    #if the component we are analysing is one of the products of the reaction, we add to the col_term and save its ID
                        col_term *= Y[comp_product.ID]/comp_product.equilibriumyield(x, mDM)
                        product.append(comp_product.ID)
            if partner.type not in SM:
                if partner.type == comp.type:
                    dY[i] += (sigma_v/2) *((Y[i]**2) - (comp.equilibriumyield(x, mDM)**2) * col_term)
                else:
                    dY[i] += - sigma_v *((Y[i] * Y[partner.ID]) - (comp.equilibriumyield(x, mDM) * partner.equilibriumyield(x, mDM)) * col_term)
                    dY[partner.ID] += - sigma_v * ((Y[i] * Y[partner.ID]) - (comp.equilibriumyield(x, mDM) * partner.equilibriumyield(x, mDM)) * col_term)
                for p in range(0, len(product)):
                    dY[(product[p])] += sigma_v * ((Y[i]*Y[partner.ID]) - (comp.equilibriumyield(x, mDM) * partner.equilibriumyield(x, mDM)) * col_term)
            if partner.type in SM:
                dY[i] += sigma_v * ((Y[i]) - comp.equilibriumyield(x, mDM) * col_term)
                for p in range(0, len(product)):
                    dY[(product[p])] += sigma_v * (Y[i] - comp.equilibriumyield(x, mDM) * col_term)
    for i in range(0, len(dY)):
    #loop to multiply all the dY elements by the 1/3H * ds/dx factor. we chose to multiply it now in an
    #attempt to simplify the numerical integration
        dY[i] = (1/(3 * H)) * dsdx * dY[i]
    print('x', x, 'dY', dY)
    return dY
                
            
        
