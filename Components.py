#module: Components
#this module sets the main characteristics for all the involved particles
#Author: Gabriela Alevi Dagrela

import math

class Component(object):
    #class to hold the main characteristics of each component
    #param label: sets a label for the component
    #param ID: sets the particle ID
    #param active: True/False, sets particle as active or not
    #param g: sets number of degrees of freedom, positive for bosons and negative for fermions
    #param mass: sets particle's rest mass
    #param decaywidth: sets particle's decay width, in its rest frame
    #param decayproducts: list with all the possible products of the particle decay
    #param br: list of all the Branching Ratios (as functions of T) of the decay reactions the particle can participate.
    #must be in the same order as the decayproducts list
    #param collisions: list of all the possible collisions the particle can participate
    #param collisionproducts: list of all the products of the possible collisions the particle can participate. must
    #be in the same order as the collisions list
    #param sigmaV: list of all the thermally averaged cross-sections for the possible collisions the particle can participate. must
    #be in the same order as the collisions and collisionproducts lists, and be a function of temperature (in GeV)

    def __init__(self, label, ID, active, g, mass, decaywidth, decayproducts, br, collision, collisionproducts, sigmaV):
        self.label = label
        self.ID = ID
        self.active = active
        self.g = g
        self.mass = mass
        self.decaywidth = decaywidth
        self.decayproducts = decayproducts
        self.br = br
        self.collision = collision
        self.collisionproducts = collisionproducts
        self.sigmaV = sigmaV

    def getBR(self, T):
        #returns the list of Branching Ratios for the particle at temperature T
        return self.br

    def getwidth(self, T):
        #returns the decay width of the particle at temperature T
        return self.decaywidth
