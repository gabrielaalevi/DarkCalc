#This module sets the model's main parameters, that define the DM and coannihilator characteristics

import numpy as np
import pyslha
from getProcesses import find_processes_col, find_decays, find_dof
from components import Component
from typing import List, Dict
import logging as logger

#Input from the user

def getModelData(paramCard : str, bsmPDGList : List[int]):


    # Get information from the param_card
    particle_data = pyslha.read(paramCard)
    massDict = particle_data['MASS']
    decaysDict = particle_data.decays
    # We need a separate method to get the quantum numbers and particle labels
    qnumbers = getQnumbers(paramCard)


    # Now restrict to the pdgs in bsmPDGList
    compDict = {}
    for pdg in bsmPDGList:
        if pdg not in qnumbers:
            logger.info(f'Particle with {pdg} not found in {paramCard}')
            continue
        particleInfo = qnumbers[pdg]
        if pdg in massDict:
            mass = massDict
        else:
            mass = 0.0
            logger.info(f'Mas for particle {pdg} not found in {paramCard}. Setting the mass to zero')
        comp = Component(label=particleInfo['label'],PDG = pdg, 
                         mass = mass, g = particleInfo['dof'])
        
        if pdg in decaysDict:
            decays = decaysDict[pdg]
            comp.setDecays(totalwidth=decays['totalwidth'], brs = decays.decays)
        
        compDict[pdg] = comp

    return compDict





def getQnumbers(paramCard : str) -> Dict[int,Dict]:

    with open(paramCard,'r') as f:
        data = f.read()

    data = data.lower()
    qnumberBlocks = data.split('block qnumbers ')[1:]
    pdgDict = {}
    for b in qnumberBlocks:
        lines = b.split('\n')
        lines = [l.strip() for l in lines[:] if l.strip()]
        lines = [l for l in lines[:] if l[0] != '#']
        header = lines[0].split('#')
        if len(header) == 2: # Found particle label
            label = header[1].strip()
            pdg = int(header[0])
        else: # Could not find label, set label to PDG
            label = header[0].strip()
            pdg = int(header[0])
        # Set default dict in case there is information missing
        pdgDict[pdg] = {'label' : label, 'spin' : 1, 'color' : 1, 
                        'anti-particle' : 0, 'dof' : 1}
        for l in lines[1:]:
            k,val = (int(x) for x in l.split()[:2])        
            if k == 1:
                pdgDict[pdg]['eCharge'] = val
            elif k == 2:
                pdgDict[pdg]['spin'] = val
            elif k == 3:
                pdgDict[pdg]['color'] = val
            elif k == 4:
                pdgDict[pdg]['anti-particle'] = val # val =0, means the particle is its anti-particle
        pdgDict[pdg]['dof'] = pdgDict[pdg]['spin']*pdgDict[pdg]['color']*(2**pdgDict[pdg]['anti-particle'])

    return pdgDict

pnames = [['DM', 'xm', 52,1],
         ['Coannihilator', 'b2', 2000005,1]] #list of the BSM particles present on the model

debug_version = False #if true, the output includes the value of the term in the Boltzmann Equation for each
#reaction, allowing for a better analysis of the result
name_file = 'model xm-b2' #name of the output file
#non-input

#list containing the SM particles
SM = ['u', 'ux', 'd', 'dx', 'c', 'cx', 's', 'sx', 't', 'tx', 'b', 'bx', 'g', 'z', 'h', 'em', 'ep', 've', 'vex', 'mup', 'mum', 'vm', 'vmx', 'tap', 'tam', 'vt', 'vtx', 'a', 'wm', 'wp']
#creating a list with the BSM particles considered in the model
BSM = list()
for i in range(0, len(pnames)):
    part = pnames[i][1]
    BSM.append(part)
#using Pyslha to read the param_card


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
        col_processes = find_processes_col(part, x, BSM, nsteps, SM)
        dof = find_dof(part, part_pdg, param_path)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], i, pnames[i][3],  part_mass, dof, x, col_processes,0, 0)
    else:
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        part = pnames[i][1]
        col_processes = find_processes_col(part, x, BSM, nsteps, SM)
        part_width = param_card.decays[part_pdg].totalwidth
        dec_reactions = find_decays(part_pdg, param_path)
        dof = find_dof(part, part_pdg, param_path)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], i, pnames[i][3], part_mass, dof, x, col_processes,part_width, dec_reactions)
    comp_names.append(comp)