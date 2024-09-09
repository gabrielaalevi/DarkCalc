import numpy as np

nsteps = 1000000
Tvalues = np.linspace(950, 9, nsteps) #interval of temperatures for the data, in GeV

#In pnames, we have a list of all the particles and some relevant information. Each particle must have a
#line in the matrix assigned to them, with 5 different elements. The elements should be, in order: label (string)
#particle type, particle PDG code (integer), ID (integer, starting from 1), active (integer, 1 if the user 
# wants the equation for said particle, 0 if the user does not want the Boltzmann equation for said particle), and in_equilibrium 
# (integer, 1 if particle begins in thermal equilibrium with the plasma, for a freeze-out model, and 
# 0 if particle's number density is virtually zero in the beggining of T, for a freeze-in model). It 
# is important to note that there should be only one particle with label 'DM', which will be used to 
# calculate x.

pnames = [['DM', 'xd', 57, 0, 1, 1],
         ['Coannihilator', 'ul', 1000002, 1, 1, 1]]

#list containing the model particles
SM = ['u', 'ux', 'd', 'dx', 'c', 'cx', 's', 'sx', 't', 'tx', 'b', 'bx', 'g', 'z', 'h', 'em', 'ep', 've', 'vex', 'mup', 'mum', 'vm', 'vmx', 'tap', 'tam', 'vt', 'vtx', 'a', 'wm', 'wp']
#creating a list with the BSM particles considered in the model
BSM = list()
for i in range(0, len(pnames)):
    part = pnames[i][1]
    BSM.append(part)