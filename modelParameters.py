import math
import numpy as np
import components

#This module sets the model's main parameters, that define the DM and mediator characteristics
Tlist = np.linspace(700, 7, 9000000) #interval of temperatures for the data, in GeV

particlesDict = [{'label': 'Mediator', #label of the first particle
                  'mass' : 700, #mass
                  'g' : 1, #degrees of freedom
                  'active': True, #if the particle is active in the model
                  'ID': 1,
                  'decaywidth': 0, #decay width for the particle, in GeV
                  'decayreactions': [{ #all the possible decay reactions for the particle, in a list of dictionaries
                      'products': ['DM', 'SM'], #the products of the decay, named using their labels, and in a list
                      'br': 0 #the branching ratio for the decay
                  },{
                      'products': ['SM', 'SM'],
                      'br': 0}],
                  'collisions': [{ #all the possible collisions for the particle, in a list of dictionaries
                      'partner' : 'DM', #partner in the collision
                      'products': ['SM', 'SM'], #the products of the collision, named using their labels, in a list
                      'sigmaV': lambda T: 10e-16 #the thermally averaged cross section for the collision
                  }, {
                      'partner': 'Mediator',
                      'products': ['DM', 'SM'],
                      'sigmaV': lambda T: 0
                  }],
                    'equilibrium': True #boolean to express if the particle is in thermal equilibrium at early times. If True, the particle is in thermal equilibrium. If false, we set its initial density as approximately zero.
                 },{
                    'label': 'DM',
                    'mass': 700.,
                    'g': 1,
                    'active': True,
                    'ID': 0,
                    'decaywidth': 0,
                    'decayreactions': 0,
                    'collisions': [{
                        'partner': 'DM',
                        'products': ['SM', 'SM'],
                        'sigmaV': lambda T: 10e-15
                    },{
                        'partner': 'Mediator',
                        'products': ['SM', 'SM'],
                        'sigmaV': lambda T: 0
                    }
                                  ],
                        'equilibrium': True
                    }]

                 