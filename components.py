#module: Components
#this module sets the main characteristics for all the involved particles

from thermal.equilibriumDensities import Neq, Yeq
from typing import List, Tuple, Dict
from scipy.interpolate import interp1d


class Component(object):
    #class to hold the main characteristics of each component

    def __init__(self, label, type, PDG, ID, mass, g):
        self.label = label
        self.type = type
        self.PDG = PDG
        self.ID = ID
        self.mass = mass
        self.g = g
        self.totalwidth = 0.0
        self.decays = None

    
    def setDecays(self, totalwidth : float, brs : Dict[List[int],float]) -> None:

        self.decays = {daughter_pids : br for daughter_pids,br in brs}
        self.totalwidth = totalwidth


    def Neq(self,T: float) -> float:
        """
        Calculates the equilibrium number density at a temperature T (in GeV) for the particle, 
        using Fermi-Dirac or Bose-Einstein Statistics.

        :param T: thermal bath temperature (in GeV)
        :param m: Particle mass
        :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

        """
        neq = Neq(T, self.mass, self.g)
        
        return neq
    
    def Yeq(self,T: float) -> float:
        """
        Calculates the equilibrium yield at a temperature T (in GeV) for the particle, 
        using Fermi-Dirac or Bose-Einstein Statistics.

        :param T: thermal bath temperature (in GeV)
        :param m: Particle mass
        :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

        """
        yeq = Yeq(T, self.mass, self.g)
        
        return yeq    
    


class CollisionProcess(object):
    """
    Holds information about a single process and its thermally average cross-section
    """


    def __init__(self,initialPDGs : List[int] = [], finalPDGs : List[int] = []):
        
        self.initialPDGs = initialPDGs
        self.finalPDGs = finalPDGs
        self.sigmaV = lambda T: 0.0

    def setSigmaV(self, Tlist : List[float], data : List[float]):

        sigmaV = interp1d(Tlist,data,fill_value=0.0,bounds_error=False)
        self.sigmaV = sigmaV



        