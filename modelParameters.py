import numpy as np
import pyslha
from sigmaV import find_processes_col, find_decays
from components import Component
from inputParameters import pnames, Tvalues

#This module sets the model's main parameters, that define the DM and coannihilator characteristics

#using Pyslha to read the param_card
param_card = pyslha.read(r'C:\Users\Gabi\Downloads\Faculdade\Diss\Code with MADDM\param_card.dat')

for i in range(0, len(pnames)):
    if pnames[i][0] == 'DM':
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        x = part_mass/Tvalues

comp_names = [] #list to hold all the components

#loop to create each component, using the information from pnames, param_card and taacs.csv
for i in range(0, len(pnames)):
    col_processes = 0
    if pnames[i][0] == 'DM':
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        part = pnames[i][1]
        col_processes = find_processes_col(part, x)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], pnames[i][3], pnames[i][4], pnames[i][5],  part_mass, 1, x, col_processes,0, 0)
    else:
        part_pdg = pnames[i][2]
        part_mass = param_card.blocks['MASS'][part_pdg]
        part = pnames[i][1]
        col_processes = find_processes_col(part, x)
        part_width = 0 #param_card.decays[part_pdg].totalwidth
        dec_reactions = 0 #find_decays(part_pdg)
        comp = Component(pnames[i][0], pnames[i][1], pnames[i][2], pnames[i][3], pnames[i][4], pnames[i][5], part_mass, 1, x, col_processes,part_width, dec_reactions)
    comp_names.append(comp)


# particlesDict = [{'label': 'Mediator', #label of the first particle
#                   'type': 'ul' #type of particle
#                   'active': True, #if the particle is active in the model
#                   'ID': 1,
#                   'decaywidth': 0, #decay width for the particle, in GeV
#                   'decayreactions': [{ #all the possible decay reactions for the particle, in a list of dictionaries
#                       'products': ['DM', 'SM'], #the products of the decay, named using their labels, and in a list
#                       'br': 0 #the branching ratio for the decay
#                   },{
#                       'products': ['SM', 'SM'],
#                       'br': 0}],
#                   'collisions': [{ #all the possible collisions for the particle, in a list of dictionaries
#                       'partner' : 'DM', #partner in the collision
#                       'products': ['SM', 'SM'], #the products of the collision, named using their labels, in a list
#                       'sigmaV': lambda T: 10e-16 #the thermally averaged cross section for the collision
#                   }, {
#                       'partner': 'Mediator',
#                       'products': ['DM', 'SM'],
#                       'sigmaV': lambda T: 0
#                   }],
#                     'equilibrium': True #boolean to express if the particle is in thermal equilibrium at early times. If True, the particle is in thermal equilibrium. If false, we set its initial density as approximately zero.
#                  },{
#                     'label': 'DM',
#                     'type': 'xd'
#                     'g': 1,
#                     'active': True,
#                     'ID': 0,
#                     'decaywidth': 0,
#                     'decayreactions': 0,
#                     'collisions': [{
#                         'partner': 'DM',
#                         'products': ['SM', 'SM'],
#                         'sigmaV': lambda T: 10e-15
#                     },{
#                         'partner': 'Mediator',
#                         'products': ['SM', 'SM'],
#                         'sigmaV': lambda T: 0
#                     }
#                                   ],
#                         'equilibrium': True
#                     }]

                 