#this module finds all the reactions a specific particle appears as the reagent, and interpolates its cross-section
import pandas as pd
import numpy as np
import pyslha

def find_processes_col(part, x, BSM, nsteps, taacs_path, processes_path):
#opening the taacs files
    processes_data = pd.read_csv(processes_path, header=None)
    taacs = pd.read_csv(taacs_path)
    collisions = list() #list to hold all the arrays of collisions for an specific particle
    sigmav_selfann = np.zeros(nsteps) #array to hold the cross section for all the self-annihilation processes
    sigmav_coann = np.zeros(nsteps) #array to hold the cross section for all co-annihilation processes with the particle right before part in the BSM list
    index_part = 0 #saves the index of the particle in the BSM list
#loop through all the processes in the taacs files
    for i in range(0, len(processes_data)):
        col_partner = '  ' #variable to hold the partner of a reaction
        products = list() #list to hold the products of each reaction
        partners_names = list() #list to hold the reagents of each reaction
        process_name = processes_data.iloc[i,1]
        reagents_name = process_name[:process_name.index('_')]
        products_name = process_name[(process_name.index('_')+1):]
        sigma_v = 0
#checking if particle participates in process i, and if it is a reagent
        if part in reagents_name:
#loop through the other particles to find the partner and the products
            for j in range(0, len(BSM)):
                part_j = BSM[j]
                if part_j == part:
                #saving the index for part in the BSM list
                    index_part = j
                if part_j in process_name:
                #checking if particle j is a reagent in reaction i. if so, saves it in the partners_list
                    if part_j in reagents_name:
                        partners_names.append(part_j)
                #checking if particle j is a product in reaction i. if so, saves it in the products list
                    if part_j in products_name:
                        products.append(part_j)
                        if products_name.count(part_j) == 2:
                        #if part_j appears twice in the products_name, we add it twice to the products list
                            products.append(part_j)
            if reagents_name.count(part) == 1:
            #if the part appears only once in the reagents_name, we find the other partner and save it in col_partner
                for k in range(0, len(partners_names)):
                    if partners_names[k] != part:
                        col_partner = partners_names[k]
            if reagents_name.count(part) == 2:
            #if the part appears twice in the reagents name, then it is colliding with itself, and we save the part as the col_partner
                col_partner = part
#reading the data from the taacs files
            x_data = taacs.iloc[:,0]
            taacs_values = taacs.iloc[:, (i+1)]
            sigma_v = np.interp(x, x_data, taacs_values)
            if not products:
#for reactions whose products are only particles from the SM
                if col_partner == part:
#self-annihilation
                    for k in range(0, nsteps):
                        sigmav_selfann[k] += sigma_v[k]
                elif col_partner == BSM[(index_part-1)]:
#coannihilation with the particle of immediate smaller index
                    for k in range(0, nsteps):
                        sigmav_coann[k] +=  sigma_v[k]
                else:
                    col = [col_partner, products, sigma_v]
                    collisions.append(col)
#adding the conversion reaction
                    col = ['SM', [col_partner], sigma_v]
                    collisions.append(col)
            else:
#if we have one or two BSM particles as products in the process
                col = [col_partner, products, sigma_v]
                collisions.append(col)
    col_selfann = [part, [], sigmav_selfann]
    collisions.append(col_selfann)
    col_coann = [BSM[(index_part-1)], [], sigmav_coann]
    collisions.append(col_coann)
    col_conversion = ['SM', [BSM[(index_part-1)], 'SM'], sigmav_coann]
    collisions.append(col_conversion)
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

def find_dof(part_name, part_pdg, param_path):
    param_card=open(param_path, 'r') #opening the param_card
    lines = param_card.readlines() #reading all the lines in the param_card
    pdg = str(part_pdg)
    for row in lines:
        string = 'BLOCK QNUMBERS ' + pdg
        if row.find(string) != -1:
            header_index = lines.index(row)
    spin_line = lines[header_index+2] #finding the index for the spin of the desired particle
    spin_string = spin_line[8] #number of total spin states 2S + 1
    spin_states = int(spin_string)
    spin = (spin_states - 1)/2 #spin for the particle
    flavour_line = lines[header_index+3] #finding the index of the line of colour rep
    flavour_string = flavour_line[8] #finding the number of flavours for the particle
    flavour_states = int(flavour_string)
    anti_line = lines[header_index+4] #finding the index of the 'particle/anti distinction' line
    anti_distinction_string = anti_line[8] #finding if the particle has an anti-particle (1) or if it is its own anti-particle (0)
    anti_distinction = int(anti_distinction_string)
    dof_total = spin_states * flavour_states
    if anti_distinction == 1:
    #if the particle has an anti-particle, we multiply the degrees of freedom by 2 to account for that
        dof_total = dof_total * 2
    if spin == 0 or spin == 1:
    #we use a convention where the degrees of freedom of bosons are negative, so we change that
        dof_total = dof_total * (-1) 
    return dof_total
    

