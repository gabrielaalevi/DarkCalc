#This module sets the model's main parameters, that define the DM and coannihilator characteristics

import numpy as np
import pyslha
from components import Component, CollisionProcess
from typing import List, Dict
import logging as logger

def getModelData(paramCard : str, bsmPDGList : List[int]) -> Dict[int,Component]:

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

def getCollisionProcesses(sigmaVfile : str) -> List[CollisionProcess]:

    
    # Get process dictionary
    with open(sigmaVfile) as f:
        comment_lines =  [l.strip() for l in f.readlines() if l.strip()[0] == '#']
        proc_lines = [l[1:].strip() for l in comment_lines if (l.count('#') == 1)]
    processDict = {}
    for l in proc_lines:    
        proc_index,proc_name,proc_pdgs = (l.split(',',2))
        proc_pdgs = proc_pdgs.replace('[','').replace(']','')
        initialPDGs,finalPDGs = proc_pdgs.split('_')
        processDict[proc_index] = {'name' : proc_name.strip(), 
                                'initialPDGs' : list(map(int, initialPDGs.split(','))), 
                                'finalPDGs' : list(map(int, finalPDGs.split(',')))}
    processes_data = np.genfromtxt(sigmaVfile, delimiter=',', skip_header=len(comment_lines), 
                                   names=True)
    collisions = []
    for proc_index,pInfo in processDict.items():
        process = CollisionProcess(initialPDGs=pInfo['initialPDGs'],
                                   finalPDGs=pInfo['finalPDGs'],
                                   name=pInfo['name'])
        process.setSigmaV(xlist=processes_data['x'],
                          sigmavList=processes_data[proc_index])

        collisions.append(process)

    return collisions

def loadInput(paramCard : str, sigmaVfile : str, dmPDG : int, bsmPDGList : List[int]):

    compDict = getModelData(paramCard,bsmPDGList=bsmPDGList)
    collisions = getCollisionProcesses(sigmaVfile)

    return compDict,collisions