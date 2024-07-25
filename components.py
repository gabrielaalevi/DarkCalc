#module: Components
#this module sets the main characteristics for all the involved particles

from auxFunc import equilibrium_yield
class Component(object):
    #class to hold the main characteristics of each component

    @classmethod
    def from_dict(label,particleDict):
        #creates object from the particleDict with the required properties

        newComponent = Component(label, active, ID, g, mass, 
                                 decayreactions, 
                                 collisions, equilibrium)
        return newComponent

    
    def __init__(self, label, active, ID, g, mass, decaywidth, decayreactions, collisions, equilibrium):
        self.label = label
        self.active = active
        self.ID = ID
        self.g = g
        self.mass = mass
        self.decaywidth = decaywidth
        self.decayreactions = decayreactions
        self.collisions = collisions
        self.equilibrium = equilibrium

    def equilibriumyield(self, x, mDM):
        Y_eq = equilibrium_yield(self.mass, x, mDM, self.g)
        return Y_eq