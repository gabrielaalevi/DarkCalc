#this module creates the boltzmann equation for each component

import math
import numpy as np
from scipy.special import kn, zeta
import AuxFunc
import Components

class BoltzmannEq(object):
    
    def Boltz(Y, x, comp_list):
        dY = np.zeros(len(comp_list)) #list with all the rhs for all the components, ordered as in the comp_list
        Y_eq = np.zeros(len(comp_list)) #list with all the equilibrium yield values for all the components, ordered as in the comp_list
        ID_list = []
        for p in range(0, len(comp_list)):
        #loop to create a list of the IDs of all particles in comp_list
            comp = comp_list[p]
            ID_list.append(comp.ID)
        comp = comp_list[0]
        T = comp.mass/x
        for k in range(0, len(comp_list)): 
        #loop to calculate the equilibrium yield of each component
            comp = comp_list[k]
            Y_eq[k] = AuxFunc.equilibrium_yield(comp.mass, T, comp.g)
        H = AuxFunc.hubblerate(T) #hubble rate at temperature T
        s = AuxFunc.entropydensity(T) #entropy density at temperature T
        dsdx = AuxFunc.dsdx(x, (comp_list[0]).mass)
        for i in range(0, len(comp_list)):
        #loop parsing through all components
            comp = comp_list[i]
            width = comp.getwidth(T)
            br = comp.getBR(T)
            mass = comp.mass
            #decay term for component i
            if comp.decayproducts != 0:
                for j in range(0, len(comp.decayproducts)):
                #loop parsing through all the possible decays of component i
                    dec_term = 1
                    product = []
                    for l in range(0, len(comp.decayproducts[j])):
                    #loop parsing through the products of decay j
                        if (comp.decayproducts[j][l]) in ID_list:
                        #this condition assures that particles that are in thermal equilibrium during the whole proccess (such as SM) will not be taken into account for this expression
                            dec_term *= Y[(comp.decayproducts[j][l])]/Y_eq[(comp.decayproducts[j][l])]
                            product.append(comp.decayproducts[j][l])
                    dY[i] += - (1/(3 * H)) * dsdx * (kn(1, mass/T)/kn(2, mass/T)) * (width(T)/s) * (Y[i] - Y_eq[i] * br[j](T) *dec_term)
                    for l in range(0, len(product)):
                        #adding the correspondent source term to the dY equation of each product
                        dY[(product[l])] += (1/(3 * H)) * dsdx * (kn(1, mass/T)/kn(2, mass/T)) * (width(T)/s) * (Y[i] - Y_eq[i] * br[j](T) * dec_term)
            #collision term for component i
            for m in range(0, len(comp.collisionproducts)):
            #loop parsing through all the collisions component i can participate in as a reagent
                col_term = 1
                product = []
                partner = comp.collision[m]
                for n in range(0, len(comp.collisionproducts[m])):
                #loop parsing through the products of the collision m
                    if (comp.collisionproducts[m][n]) in ID_list:
                        col_term *= Y[comp.collisionproducts[m][n]]/Y_eq[(comp.collisionproducts[m][n])]
                        product.append(comp.collisionproducts[m][n])
                if partner in ID_list:
                    dY[i] += - (1/(3 * H)) * dsdx* comp.sigmaV[m](T) *((Y[i] * Y[partner])/(Y_eq[i] * Y_eq[partner]) - col_term)
                    dY[partner] += - (1/(3*H)) * dsdx * comp.sigmaV[m](T) * ((Y[i] * Y[partner])/(Y_eq[i] * Y_eq[partner]) - col_term)
                    for p in range(0, len(product)):
                        dY[(product[p])] += (1/(3 * H)) * dsdx * comp.sigmaV[m](T) * ((Y[i] * Y[partner])/(Y_eq[i] * Y_eq[partner]) - col_term)
                if partner not in ID_list:
                    dY[i] += (1/(3*H)) * dsdx * comp.sigmaV[m](T) * ((Y[i]/ Y_eq[i]) - col_term)
                    for p in range(0, len(product)):
                        dY[(product[p])] += (1/(3 * H)) * dsdx * comp.sigmaV[m](T) * ((Y[i]/Y_eq[i]) - col_term)
        print (dY)
        return dY
                    
                
            
    