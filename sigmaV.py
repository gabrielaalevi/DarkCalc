#this module finds all the reactions a specific particle appears as the reagent, and interpolates its cross-section
import pandas as pd
import numpy as np
import pyslha

def find_processes_col(part, x, BSM, nsteps):
    processes_data = pd.read_csv("processes_taacs.csv", header=None)
    taacs = pd.read_csv("taacs.csv")
    collisions = list() #list to hold all the arrays of collisions for an specific particle
    sigmav_selfann = np.zeros(nsteps) #array to hold the cross section for all the self-annihilation processes, with the intention of simplyfing the integration
    sigmav_coann = np.zeros(nsteps) #array to hold the cross section for all co-annihilation processes with the particle right before part in the BSM list
    index_part = 0 #saves the index of part in the BSM list
#loop through all the processes in the taacs files
    for i in range(0, len(processes_data)):
        col_partner = '  ' #variable to hold the partner of a reaction
        products = list() #list to hold the products of each reaction
        partners_names = list() #list to hold the reagents of each reaction
        process_name = processes_data.iloc[i,1]
        sigma_v = 0
#checking if particle participates in process i
        if part in process_name and process_name.index(part) < process_name.index('_'):
#loop through the other particles to find the partner and the products
            for j in range(0, len(BSM)):
                part_j = BSM[j]
                if part_j == part:
                    index_part = j
                if part_j in process_name:
#checking if particle j is a reagent in reaction i. if so, saves it in the partners_list
                    if process_name.index(part_j) < process_name.index('_'):
                        partners_names.append(part_j)
#checking if particle j is a product in reaction i. if so, saves it in the products list
                    if process_name.index(part_j) > process_name.index('_'):
                        products.append(part_j)
#checking if the partner is the particle analysed (self-annihilation) or not (co-annihilation)
            if process_name.count(part) == 1:
                for k in range(0, len(partners_names)):
                    if partners_names[k] != part:
                        col_partner = partners_names[k]
            if process_name.count(part) == 2:
                col_partner = part
            x_data = taacs.iloc[:,0]
            taacs_values = taacs.iloc[:, (i+1)]
            sigma_v = np.interp(x, x_data, taacs_values)
            if not products:
                if col_partner == part:
                    for k in range(0, nsteps):
                        sigmav_selfann[k] += sigma_v[k]
                elif col_partner == BSM[(index_part-1)]:
                    for k in range(0, nsteps):
                        sigmav_coann[k] += sigma_v[k]
                else:
                    col = [col_partner, products, sigma_v]
                    collisions.append(col)
            else:
                col = [col_partner, products, sigma_v]
                collisions.append(col)
    col_selfann = [part, [], sigmav_selfann]
    collisions.append(col_selfann)
    col_coann = [BSM[(index_part-1)], [], sigmav_coann]
    collisions.append(col_coann)
    print(collisions)
    return collisions

def find_decays(part, param_path):
#function to find all the decay reactions for a particle
    decayreactions = list() #list to hold all the decay reactions
    params = pyslha.read(param_path) #reading the param_card
    decay = params.decays[part].decays
    for i in range(0, len(decay)):
    #parsing through all the decays for part
        dec = decay[i]
        daughters = dec.ids
        dec_reaction = [daughters, dec.br]
        decayreactions.append(dec_reaction)
    print(decayreactions)
    return decayreactions
        

