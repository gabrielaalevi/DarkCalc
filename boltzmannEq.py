#this module creates the boltzmann equation for each component

import numpy as np
from scipy.special import kn
import auxFunc

def boltz(x, Y, comp_names, SM, mDM, x_new):
    dY = np.zeros(len(comp_names)) #list with all the rhs for all the components, ordered as in the comp_list
    H = auxFunc.hubblerate(x, mDM) #hubble rate at temperature T
    s = auxFunc.entropydensity(x, mDM) #entropy density at temperature T
    dsdx = auxFunc.dsdx(x, mDM) #variation of entropy with x
    T = mDM/x
    for i in range(0, len(comp_names)):
    #loop parsing through all components
        comp = comp_names[i]
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
                dY[i] += (1/(3 * H)) * dsdx * (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM) *dec_term)
                a = (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM) *dec_term)
                print('comp', comp.type, 'products', product, 'a', a)
                for l in range(0, len(product)):
                #adding the correspondent source term to the dY equation of each product
                    dY[(product[l])] += - (1/(3 * H)) * dsdx * (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[i] - comp.equilibriumyield(x, mDM)* dec_term)
        
    #collision term for component i
        for m in range(0, len(comp.collisions)):
        #loop parsing through all the collisions component i can participate in as a reagent
            col_term = 1
            prod_yeqnonrelativ = 1
            product = []
            partner = comp.collisions[m][0] #saving the TYPE of the partner
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
                        col_term *= Y[comp_product.ID]/comp_product.equilibriumyield(x, mDM)
                        prod_yeqnonrelativ = prod_yeqnonrelativ * comp_product.yeqnonrelativ(x, mDM)
                        product.append(comp_product.ID)
            if partner != 'SM':
            #for self-annihilation, co-annihilation, and double conversion interactions
                sigma_v = ((comp.collisions[m][2][index])* (comp.yeqnonrelativ(x, mDM) * partner_comp.yeqnonrelativ(x, mDM)))/((comp.equilibriumyield(x, mDM) * partner_comp.equilibriumyield(x, mDM)))
                dY[i] += (1/(3 * H)) * dsdx * sigma_v * ((Y[i] * Y[partner_comp.ID]) - (comp.equilibriumyield(x, mDM) * partner_comp.equilibriumyield(x, mDM) * col_term))
                a = (1/(3 * H)) * dsdx
                b = (1/(3 * H)) * dsdx * sigma_v * ((Y[i] * Y[partner_comp.ID]) - (comp.equilibriumyield(x, mDM) * partner_comp.equilibriumyield(x, mDM) * col_term))
                print('x', x, 'comp', comp.type, 'partner', partner, 'sigv', sigma_v, 'Y_comp', Y[i], 'Y_eq comp', comp.equilibriumyield(x, mDM), 'Y_part', Y[partner_comp.ID], 'Y_eq part', partner_comp.equilibriumyield(x, mDM), 'col term', col_term, '1/3h dsdx', a)
                print(b)
            else:
            #for conversion interaction
                sigma_v = ((comp.collisions[m][2][index])* (comp.yeqnonrelativ(x, mDM) *prod_yeqnonrelativ))/ (comp.equilibriumyield(x, mDM))
                dY[i] += (1/(3 * H)) * dsdx * sigma_v * ((Y[i]) - comp.equilibriumyield(x, mDM) * col_term)
                a = (1/(3 * H)) * dsdx
                b = (1/(3 * H)) * dsdx * sigma_v * ((Y[i]) - comp.equilibriumyield(x, mDM) * col_term)
                print('x', x, 'comp', comp.type, 'partner', partner, 'sigv', sigma_v, 'Y_comp', Y[i], 'Y_eq comp', comp.equilibriumyield(x, mDM), 'Y_part', Y[partner_comp.ID], 'Y_eq part', partner_comp.equilibriumyield(x, mDM), 'col term', col_term, '1/3h dsdx', a)
                print(b)
                for p in range(0, len(product)):
                    dY[product[p]] += -(1/(3 * H)) * dsdx * sigma_v * ((Y[i]) - comp.equilibriumyield(x, mDM) * col_term)
    print('x', x, 'dY', dY)
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
    H = auxFunc.hubblerate(x, mDM) #hubble rate at temperature T
    dsdx = auxFunc.dsdx(x, mDM) #variation of entropy with x
    s = auxFunc.entropydensity(x, mDM) #entropy density at temperature T
    tol = abs(x - x_new[0])
    T = mDM/x
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
                reaction_terms[counter] = ((1/(3 * H))* dsdx * (kn(1, comp.mass/T)/kn(2, comp.mass/T)) * (comp.decaywidth/s) * (Y[2*i] - Y[(2*i+1)] *dec_term))/H
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
                
            
        
