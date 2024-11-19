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
    

class CollisionProcesses(object):


    def __init__(self):
        
        self._processes = {}
    
    def addProcess(self, Tlist : List[float], collisionData : List[float], 
                    initialStates : List["Component"], finalStates : List["Component"]) -> None:
        
        initialPDGs = sorted([comp.PDG for comp in initialStates])
        finalPDGs = sorted([comp.PDG for comp in finalStates])
        sigmaF = interp1d(Tlist,collisionData,fill_value=0.0,bounds_error=False)
        if not initialPDGs in self._processes:
            self._processes[initialPDGs] = {}
        self._processes[initialPDGs][finalPDGs] = sigmaF

    def sigmaV(self,T: float, 
                initialStates : List["Component"], 
                finalStates : List["Component"]) -> float:

        initialPDGs = sorted([comp.PDG for comp in initialStates])
        if not initialPDGs in self._processes:
            return None
        finalPDGs = sorted([comp.PDG for comp in finalStates])
        if not finalPDGs in self._processes[initialPDGs]:
            return None
        
        sv = self._processes[initialPDGs][finalPDGs](T)

        return sv
        

        