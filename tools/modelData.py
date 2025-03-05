#module: Components
#this module sets the main characteristics for all the involved particles

from thermal.equilibriumDensities import Neq, Yeq, gstar
from typing import List, Dict, Optional
from scipy.interpolate import interp1d
from scipy.special import kn
from tools.logger import logger
import numpy as np
import pyslha
from typing import List, Dict, Union
from tools.logger import logger
import os




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

    def __str__(self) -> str:

        if self.label:
            return self.label
        else:
            return str(self.PDG)
        
    def __repr__(self) -> str:
        return str(self)
        

    
    def setDecays(self, totalwidth : float, decays : list, pdgList : List[int]) -> None:

        
        self.decays = {}
        for decay in decays:
            finalIDs = decay.ids
            finalIDs = np.abs(finalIDs).tolist()
            finalIDs = [fid if fid in pdgList else 0 for fid in finalIDs]
            finalIDs = sorted(finalIDs)
            self.decays[tuple(finalIDs)] = decay.br
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
   
    def gammaInv(self,T: float) -> float:
        """
        Computes the inverse of the relativistic factor (1/gamma), which is used
        when computing effective decay rates.
        For small T or large T (compared to the component's mass) it is useful to expand
        the factor to avoid numerical instabilities.
        """

        x = self.mass/T
        if x > 10: # relativistic regime
            return 1.0 - 3.0/(2*x) + 15./(8*x**2)
        elif x < 0.05: # non-relativistic regime
            return x/2.0
        else:
            return kn(1,x)/kn(2,x)


class CollisionProcess(object):
    """
    Holds information about a single process and its thermally average cross-section
    """


    def __init__(self,initialPDGs : List[int] = [], finalPDGs : List[int] = [], 
                 name : Optional[str] = "") -> None:
        """
        Stores the thermally averaged cross-section for the process initialPDGs -> finalPDGs.
        Only the absolute values of PDGs are used, so any charge-conjugated process will
        result in the same CollisionProcess.
        """
        
        self.initialPDGs = sorted(np.abs(initialPDGs).tolist())
        self.finalPDGs = sorted(np.abs(finalPDGs).tolist())
        self.sigmaV = lambda x: 0.0
        self.name = name

    def __str__(self):
        
        proc_str = f'{self.name}:{tuple(self.initialPDGs)} -> {tuple(self.finalPDGs)}'

        return proc_str
    
    def __repr__(self):
        return str(self)

    def changesPDG(self,pdg : int) -> bool:
        """
        Checks if the process changes the Yield for particle.
        If the particle's pdg does not appear on the process or
        if the number of particles in the initial state equals the number
        in the final state, return False. Return True otherwise.

        :param pdg: Particle's pdg code.
        """

        if pdg not in (self.initialPDGs+self.finalPDGs):
            return False
        elif self.initialPDGs.count(pdg) == self.finalPDGs.count(pdg):
            return False
        else:
            return True

    def setSigmaV(self, xlist : List[float], sigmavList : List[float]):

        sigmaV = interp1d(xlist,sigmavList,fill_value="extrapolate",bounds_error=False)
        self.sigmaV = sigmaV



class ModelData(object):
    """
    Holds all the necessary information about the model, its components (BSM particles)
    and collision processes
    """

    def __init__(self, dmPDG : Union[int,str],  
                 bsmPDGList : Union[List[int],List[str]], 
                 inputFile : Optional[str] = None):
        """
        :param bsmPDGList: Used to selected the BSM particles to include in the Boltzmann equations. 
                            Currently only a Dark Matter and a Mediator are allowed.
        :param dmPDG: Used to define the PDG particle
        """
        
        if dmPDG not in bsmPDGList:
            bsmPDGList.append(dmPDG)

        if len(bsmPDGList) > 2:
            logger.error(f"The Model can only contain up to two BSM Particles (not {len(bsmPDGList)}!)")
            raise ValueError()
        elif 0 in bsmPDGList:
            logger.error(f"The Model can not contain BSM Particles with PDG = 0!")
            raise ValueError()
        # elif any(pdg < 0 for pdg in bsmPDGList):
        #     logger.error(f"The BSM particles should be defined with positive PDGs!")
        #     raise ValueError()
        
        smComponent = Component(label='SM',PDG=0, mass = 0.0, g = 1.0, ID = 0)
        # The a dynamical number of degrees of freedom for the SM
        smComponent.g = gstar
        self.componentsDict : Dict[int,Component] = {0 : smComponent}
        self.collisionProcesses : List[CollisionProcess] = []
        self.dmPDG = dmPDG # Dark Matter PDG code
        self.pdgList = bsmPDGList[:]

        if inputFile is not None:
            self.getParticleDataFrom(inputFile)
            self.getCollisionProcessesFrom(inputFile)

    def __str__(self) -> str:

        dmLabel = self.componentsDict[self.dmPDG].label
        smLabel = self.componentsDict[0].label
        otherBSMLabels = [comp.label for comp in self.componentsDict.values() 
                            if comp.label not in [dmLabel,smLabel]]
        mStr = f'Model : {smLabel} + {dmLabel}(DM)'
        if otherBSMLabels:
            mStr += ' + ' + ' + '.join(otherBSMLabels)

        return mStr
    
    def __repr__(self) -> str:
        return str(self)
    
    @classmethod
    def loadModel(cls, parser : dict, outputFile: str) -> "ModelData":

        logger.debug(f'Loading model from {outputFile}')
        if not os.path.isdir(outputFile):
            logger.error(f'Output file {outputFile} not found')
            return False
        

        dm = parser['Model']['darkmatter']
        if 'bsmParticles' in parser['Model']:
            bsmList = str(parser['Model']['bsmParticles']).split(',')
        else:
            bsmList = []
        model = ModelData(dmPDG=dm, bsmPDGList=bsmList, 
                          inputFile=outputFile)
        logger.info(f'Successfully loaded {model}')

        return model

    def getParticleDataFrom(self, inputFile : str):

        # Get information from the param_card

        with open(inputFile,'r') as f:
            output_data = f.read()
        paramCard_data = output_data.split('<slha>')[1].split('</slha>')[0]

        particle_data = pyslha.read(paramCard_data)
        massDict = particle_data.blocks['MASS']
        decaysDict = particle_data.decays
        # We need a separate method to get the quantum numbers and particle labels
        qnumbers = self.getQnumbersFrom(paramCard)
        labels2pdgs = {qnumbers[pdg]['label'] : pdg for pdg in qnumbers}

        # Convert labels to PDGs if labels were given as input
        self.pdgList = [labels2pdgs[pdg] if pdg in labels2pdgs else int(pdg)
                        for pdg in self.pdgList[:]]
        if self.dmPDG in labels2pdgs:
            self.dmPDG = labels2pdgs[self.dmPDG]
        else:
            self.dmPDG = int(self.dmPDG)

        # Now restrict to the pdgs in bsmPDGList
        for pdg in self.pdgList:
            if pdg not in qnumbers:
                logger.info(f'Particle with {pdg} not found in {paramCard}')
                continue
            particleInfo = qnumbers[pdg]
            if pdg in massDict:
                mass = massDict[pdg]
            else:
                mass = 0.0
                logger.info(f'Mas for particle {pdg} not found in {paramCard}. Setting the mass to zero')
            comp = Component(label=particleInfo['label'],PDG = pdg, 
                            mass = mass, g = particleInfo['dof'],
                            ID = len(self.componentsDict))
            if pdg in decaysDict:
                decays = decaysDict[pdg]
                comp.setDecays(totalwidth=decays.totalwidth, decays = decays.decays, 
                               pdgList = self.pdgList)
            
            self.componentsDict[pdg] = comp
        
        if self.dmPDG not in self.componentsDict:
            logger.error(f"Dark Matter (PDG = {self.dmPDG}) not found in {paramCard}")
            raise ValueError()
        
        # Set Dark Matter mass to use for scaling the temperature
        self.mDM = self.componentsDict[self.dmPDG].mass

        self.createLabel2PDGDict()

        # If collisions have been already loaded, make sure
        # that all the PDGs not appearing in self.pdgList are in thermal equilibrium (SM)
        # and set their PDG to zero
        if self.collisionProcesses:
           for proc in self.collisionProcesses:
               initPDGs = proc.initialPDGs[:]
               proc.initialPDGs = [x if x in self.pdgList else 0 for x in initPDGs ]
               finalPDGs = proc.finalPDGs[:]
               proc.finalPDGs = [x if x in self.pdgList else 0 for x in finalPDGs ]

    def getQnumbersFrom(self,paramCard : str) -> Dict[int,Dict]:

        with open(paramCard,'r') as f:
            data = f.read()

        data = data.lower()
        qnumberBlocks = data.split('block qnumbers ')[1:]
        pdgDict = {}
        for b in qnumberBlocks:
            b = b[:b.find('\ndecay')]
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
            if pdgDict[pdg]['spin']%2 == 0:
                pdgDict[pdg]['dof'] = -pdgDict[pdg]['dof'] # Set g < 0 for fermions

        return pdgDict

    def getCollisionProcessesFrom(self,sigmaVfile : str):
        
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
        
        for proc_index,pInfo in processDict.items():
            process = CollisionProcess(initialPDGs=pInfo['initialPDGs'],
                                    finalPDGs=pInfo['finalPDGs'],
                                    name=pInfo['name'])
            process.setSigmaV(xlist=processes_data['x'],
                            sigmavList=processes_data[proc_index])

            self.collisionProcesses.append(process)

        # Assume all the PDGs not appearing in self.pdgList are in thermal equilibrium (SM)
        # and set their PDG to zero
        if self.pdgList:
           for proc in self.collisionProcesses:
               initPDGs = proc.initialPDGs[:]
               proc.initialPDGs = [x if x in self.pdgList else 0 for x in initPDGs ]
               finalPDGs = proc.finalPDGs[:]
               proc.finalPDGs = [x if x in self.pdgList else 0 for x in finalPDGs ]

    def createLabel2PDGDict(self) -> None:
        if not hasattr(self,'_label2pdgDict'):
            self._label2pdgDict = {comp.label : comp.PDG 
                         for comp in self.componentsDict.values()}
            self._label2pdgDict.update({comp.PDG : comp.PDG 
                         for comp in self.componentsDict.values()})
   
    def convert2PDG(self,label : str) -> int:
       
        return self._label2pdgDict[label]
    
