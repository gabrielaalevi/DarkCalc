#This module sets the model's main parameters, that define the DM and coannihilator characteristics

import numpy as np
import pyslha
from getProcesses import find_processes_col, find_decays, find_dof
from components import Component

#Input from the user

nsteps = 1000000 #number of values of x to be used
T_i = 555 #initial temperature, in GeV
T_f = 5 #final temperature, in GeV
taacs_path = r'C:\Users\Gabi\Downloads\Faculdade\Diss\Code with MADDM\taacs.csv' #path to the taacs.csv
processes_path = r'C:\Users\Gabi\Downloads\Faculdade\Diss\Code with MADDM\processes_taacs.csv' #path to the processes_taacs.csv
param_path = r'C:\Users\Gabi\Downloads\Faculdade\Diss\Code with MADDM\param_card.dat' #path to the param_card

#In pnames, we have a list of all the particles and some relevant information. Each particle must have a
#line in the matrix assigned to them, with 5 different elements. The elements should be, in order: label (string)
#particle type, particle PDG code (integer), in_equilibrium (integer, 1 if particle begins in thermal equilibrium 
# with the plasma, for a freeze-out model, and 0 if particle's number density is virtually zero in the 
# beggining of T, for a freeze-in model). It is important to note that there should be only one particle with 
# label 'DM', which will be used to calculate the value of x.

pnames = [['DM', 'xm', 52,1],
         ['Partner', 'b2', 2000005,1]] #list of the BSM particles present on the model

debug_version = True #if true, the output includes the value of each term in the Boltzmann Equation for each
#reaction, allowing for a better analysis of the result
name_file = 'model xm-b2' #name of the output file



#non-input

Tvalues = np.linspace(T_i, T_f, nsteps) #interval of temperatures for the data, in GeV
#list containing the SM particles
SM = ['u', 'ux', 'd', 'dx', 'c', 'cx', 's', 'sx', 't', 'tx', 'b', 'bx', 'g', 'z', 'h', 'em', 'ep', 've', 'vex', 'mup', 'mum', 'vm', 'vmx', 'tap', 'tam', 'vt', 'vtx', 'a', 'wm', 'wp']
#creating a list with the BSM particles considered in the model
BSM = list()
for i in range(0, len(pnames)):
    part = pnames[i][1]
    BSM.append(part)
#using Pyslha to read the param_card
param_card = pyslha.read(param_path)

#creating the array of x, using the DM mass and the Tvalues
for i in range(0, len(pnames)):
    if pnames[i][0] == 'DM':
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        x = part_mass/Tvalues
        mDM = part_mass

comp_names = [] #list to hold all the components
#loop to create each component, using the information from pnames, param_card and taacs.csv
for i in range(0, len(pnames)):
    col_processes = 0
    if pnames[i][0] == 'DM':
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        part = pnames[i][1]
        col_processes = find_processes_col(part, x, BSM, nsteps, taacs_path, processes_path)
        dof = find_dof(part, part_pdg, param_path)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], i, pnames[i][3],  part_mass, dof, x, col_processes,0, 0)
    else:
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        part = pnames[i][1]
        col_processes = find_processes_col(part, x, BSM, nsteps, taacs_path, processes_path)
        part_width = param_card.decays[part_pdg].totalwidth
        dec_reactions = find_decays(part_pdg, param_path)
        print(dec_reactions)
        dof = find_dof(part, part_pdg, param_path)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], i, pnames[i][3], part_mass, dof, x, col_processes,part_width, dec_reactions)
    comp_names.append(comp)