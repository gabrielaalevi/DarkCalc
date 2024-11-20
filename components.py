#module: Components
#this module sets the main characteristics for all the involved particles

from thermal.equilibriumDensities import Neq, Yeq
from typing import List, Dict, Optional
from scipy.interpolate import interp1d


class Component(object):
    #class to hold the main characteristics of each component

    def __init__(self, label, PDG, ID, mass, g):
        self.label = label
        self.PDG = PDG
        self.ID = ID
        self.mass = mass
        self.g = lambda T : g # For convenience set it as a dummy function
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
        neq = Neq(T, self.mass, self.g(T))
        
        return neq
    
    def Yeq(self,T: float) -> float:
        """
        Calculates the equilibrium yield at a temperature T (in GeV) for the particle, 
        using Fermi-Dirac or Bose-Einstein Statistics.

        :param T: thermal bath temperature (in GeV)
        :param m: Particle mass
        :param g: Particle's degrees of freedom. Use g < 0 for bosons and g > 0 for fermions

        """
        yeq = Yeq(T, self.mass, self.g(T))
        
        return yeq    
    


class CollisionProcess(object):
    """
    Holds information about a single process and its thermally average cross-section
    """


    def __init__(self,initialPDGs : List[int] = [], finalPDGs : List[int] = [], 
                 name : Optional[str] = "") -> None:
        
        self.initialPDGs = sorted(initialPDGs)
        self.finalPDGs = sorted(finalPDGs)
        self.sigmaV = lambda x: 0.0
        self.name = name

    def setSigmaV(self, xlist : List[float], sigmavList : List[float]):

        sigmaV = interp1d(xlist,sigmavList,fill_value=0.0,bounds_error=False)
        self.sigmaV = sigmaV

    def __str__(self):
        
        proc_str = f'{self.name}:{tuple(self.initialPDGs)} -> {tuple(self.finalPDGs)}'

        return proc_str
    
    def __repr__(self):
        return str(self)



        