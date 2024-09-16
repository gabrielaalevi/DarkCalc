#module: Components
#this module sets the main characteristics for all the involved particles

from auxFunc import equilibrium_yield
class Component(object):
    #class to hold the main characteristics of each component

    def __init__(self, label, type, PDG, ID, in_equilibrium, mass, g, xarray, collisions, decaywidth, decayreactions):
        self.label = label
        self.type = type
        self.PDG = PDG
        self.ID = ID
        self.in_equilibrium = in_equilibrium
        self.mass = mass
        self.g = g
        self.xarray = xarray
        self.collisions = collisions
        self.decaywidth = decaywidth
        self.decayreactions = decayreactions

    def equilibriumyield(self, x, mDM):
        Y_eq = equilibrium_yield(self.mass, x, mDM, self.g)
        return Y_eq