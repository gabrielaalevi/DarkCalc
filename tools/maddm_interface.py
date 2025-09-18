import os,shutil
from tools.logger import logger
import tempfile
import subprocess
import numpy as np
import json
from io import StringIO
from tools.modelData import smPDGs
from typing import Union


def runMadDM(parser : dict) -> str:
    """
    Run MadDM to compute widths and relevant collision cross-sections.
    
    :param parser: Dictionary with parser sections.
    
    :return: True if successful. Otherwise False.
    """
        
    #Get run folder:    
    outputFolder = os.path.abspath(parser['Options']['outputFolder'])
    maddmFolder = os.path.join(outputFolder,'maddm')
    if not os.path.isdir(maddmFolder):
        os.makedirs(maddmFolder)    
    
    modelDir = os.path.abspath(parser['Model']['modelDir'])
    if not os.path.isdir(modelDir):
        # Check for model restriction in the name
        mDir = modelDir.rsplit('-',1)[0]
        if not os.path.isdir(mDir):
            logger.error(f'Model folder {modelDir} (or {mDir}) not found')
            raise ValueError()

    dm = parser['Model']['darkmatter']
    if 'bsmParticles' in parser['Model']:
        bsmList = str(parser['Model']['bsmParticles']).split(',')
        # Make sure the dmPDG does not appear twice
        bsmList = [p for p in bsmList[:] if p != dm]
    else:
        bsmList = []

    addConversion = parser['Options'].get('addConversion',False)
    # If there are no additional BSM particles, addConversion should be False
    if not bsmList:
        addConversion = False

    #Generate commands file:       
    commandsFile,cFilePath = tempfile.mkstemp(suffix='.txt', 
                                              prefix='maddm_commands_', 
                                              dir=outputFolder)    
    os.close(commandsFile)
    commandsFileF = open(cFilePath,'w')
    commandsFileF.write(f'import model {modelDir}\n')
    commandsFileF.write(f'define darkmatter {dm}\n')
    for p in bsmList:
        commandsFileF.write(f'define coannihilator {p}\n')
    commandsFileF.write('generate relic_density\n')
    commandsFileF.write(f'output {maddmFolder}\n')
    commandsFileF.write('launch\n')
    if 'paramCard' in parser['Model']:
        paramCard = os.path.abspath(parser['Model']['paramCard'])
        commandsFileF.write(f'{paramCard} \n')
    comms = parser["SetParameters"]
    #Set model parameters
    for key,val in comms.items():
        commandsFileF.write(f'set {key} {val}\n')
    
    
    if 'computeWidths' in parser['Model']:
        pList = str(parser['Model']['computeWidths']).split(',')
        pStr = ' '.join(pList)
        commandsFileF.write(f'compute_widths {pStr}\n')
        commandsFileF.write('done\n')
    commandsFileF.close()
    mg5Folder = parser['Options']['MadGraphPath']
    mg5Folder = os.path.abspath(mg5Folder)        
    if not os.path.isfile(os.path.join(mg5Folder,'bin','maddm.py')):
        logger.error(f'Executable maddm.py not found in {mg5Folder}')
        raise FileNotFoundError()
    # Comput widths
    with open(cFilePath, 'r') as f: 
        logger.debug(f'Running MadDM with commands:\n {f.read()} \n')
    run = subprocess.Popen(f'./bin/maddm.py -f {cFilePath}',shell=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                cwd=mg5Folder)
        
        
         
    output,errorMsg = run.communicate()
    logger.debug(f'MadDM process error:\n {errorMsg.decode("utf-8")} \n')
    logger.debug(f'MadDM process output:\n {output.decode("utf-8")} \n')

    if os.path.isfile(cFilePath):
        os.remove(cFilePath)

    # Check if cross-sections were generated and saved to taacs.csv
    sigmaVFile = os.path.join(maddmFolder,'output','taacs.csv')
    if not os.path.isfile(sigmaVFile):
        logger.error(f"Error computing sigmaV with MadDM ({sigmaVFile} not found)")
        raise FileNotFoundError()
    else:
        return mergeOutput(outputFolder, maddmFolder, addConversion)

@np.vectorize
def is_sm(pdg : int) -> bool:
    """
    Checks if the PDG belongs to the SM
    """
    return abs(pdg) in smPDGs


def simplifyProcess(pDict : dict) -> dict:
    """
    Simplify the collision process. The SM PDGs are replaced by zero
    and if the process is of the type BSM_i BSM_j <-> SM_a SM_b,
    its name is simplified to BSM_i BSM_j <-> SM SM.
    It also replaces all PDGs by their absolute values.
    """

    pDict_new = {k : v for k,v in pDict.items()}
    pname = pDict['name']
    initialPDGs = pDict['initialPDGs']
    finalPDGs = pDict['finalPDGs']
    if sum(is_sm(initialPDGs)) == len(initialPDGs):
        in_new = 'SMSM'
    else:
        in_new = pname.split('_')[0]
    if sum(is_sm(finalPDGs)) == len(finalPDGs):
        out_new = 'SMSM'
    else:
        out_new = pname.split('_')[1]
    pname_new = in_new+'_'+out_new
    initialPDGs_new = np.where(is_sm(initialPDGs),0,initialPDGs)
    finalPDGs_new = np.where(is_sm(finalPDGs),0,finalPDGs)
    
    pDict_new['name'] = pname_new
    pDict_new['initialPDGs'] = sorted([abs(pdg) for pdg in initialPDGs_new])
    pDict_new['finalPDGs'] = sorted([abs(pdg) for pdg in finalPDGs_new])

    return pDict_new


def convertProcess(pDict : dict, index_str : str = "") -> Union[None,dict]:
    """
    Convert the collision process to a conversion process.
    It is only applicable to processes of type BSM_i BSM_j <-> SM SM,
    otherwise it returns None.
    It assumes the thermally averaged cross-section for BSM_i BSM_j <-> SM SM
    is the same for SM BSM_j <-> SM BSM_i.    
    """
    inPDGs = sorted(pDict['initialPDGs'])
    outPDGs = sorted(pDict['finalPDGs'])
    # First check if the process is of type  BSM_i BSM_j <-> SM SM
    if all(is_sm(pdg) for pdg in inPDGs) and all(not is_sm(pdg) for pdg in outPDGs):
        # Convert SM SM -> BSM_i BSM_j to SM BSM_i -> BSM_j SM
        inPDGs_new = inPDGs[:-1] + [outPDGs[0]]
        outPDGs_new = outPDGs[1:] + [inPDGs[-1]]        
    elif all(is_sm(pdg) for pdg in outPDGs) and all(not is_sm(pdg) for pdg in inPDGs):
        # Convert BSM_i BSM_j <-> SM SM to SM BSM_j <-> BSM_i SM
        inPDGs_new = inPDGs[1:] + [outPDGs[0]]
        outPDGs_new = outPDGs[1:] + [inPDGs[0]]
    else:
        return None
    
    pDict_new = {k : v for k,v in pDict.items()}
    pDict_new['initialPDGs'] = sorted(inPDGs_new)
    pDict_new['finalPDGs'] = sorted(outPDGs_new)
    pDict_new['name'] = 'conversion'
    if index_str:
        pDict_new['index'] = index_str

    return pDict_new

def mergeOutput(outputFolder : str, maddmFolder : str, addConversion: bool) -> str:
    """
    Combines the param_card and the sigmaVFile into a single file,
    similar to the MadGraph banner.
    Returns the new file name
    """

    paramCard = os.path.join(maddmFolder,'Cards','param_card.dat')
    if  not os.path.isfile(paramCard):
        logger.error(f"Param card ({paramCard} not found)")
        raise FileNotFoundError()
    sigmaVFile = os.path.join(maddmFolder,'output','taacs.csv')
    if not os.path.isfile(sigmaVFile):
        logger.error(f"SigmaV file ({sigmaVFile} not found)")
        raise FileNotFoundError()
    
    # Loads the MadDM output (taacs.csv)
    with open(sigmaVFile, 'r') as f:
        lines = f.readlines()
    data = [l for l in lines if l.strip()[0] != '#']
    processInfo = [l for l in lines if (l.strip()[0] == '#') and (l.count('#') == 1)]
    data = np.genfromtxt(data,delimiter=',',comments='#',names=True)

    # Extract the cross-sections for each process and simplify the processes
    processDict = {}
    for line in processInfo:    
        proc_index, proc_name, proc_pdgs = line.replace('\n','').split(',',2)
        proc_index = proc_index.replace('#', '').strip()
        proc_name = proc_name.strip()
        initialPDGs,finalPDGs = proc_pdgs.split('_')
        initialPDGs = sorted(json.loads(initialPDGs)) # Convert to list
        finalPDGs = sorted(json.loads(finalPDGs)) # Convert to list
        if proc_index not in data.dtype.names:
            loggger.error(f'Process {proc_index} ({proc_name}) not found in data columns.')
            raise ValueError()
        pDict = {'index' : proc_index,
                'name' : proc_name,
                'initialPDGs' : initialPDGs, 
                'finalPDGs' : finalPDGs,
                'data' : np.array(data[['x',proc_index]].tolist())
                }
        # Simplify the process (replace SM PDGs by zero and simplify process name)
        pDict = simplifyProcess(pDict)
        proc = (tuple(pDict['initialPDGs']),tuple(pDict['finalPDGs']))
        # If the process has already appeared, add its cross-section:
        if proc in processDict:
            processDict[proc]['data'][:,1] +=  pDict['data'][:,1]
        else:
            processDict[proc] = simplifyProcess(pDict)
            
    # If required, use processes of type BSM_i BSM_j <-> SM_a SM_b
    # to add conversion process (SM BSM_j <-> SM BSM_i)
    max_index = max([int(pDict['index']) for pDict in processDict.values()])
    max_index += 1
    if addConversion:
        for proc,pDict in list(processDict.items()):
            pDict_new = convertProcess(pDict,index_str=str(max_index))
            if pDict_new is None:
                continue
            proc_new = (tuple(pDict_new['initialPDGs']),
                        tuple(pDict_new['finalPDGs']))
            # Conversion process was already present, skip
            if proc_new in processDict:
                continue
            else:
                processDict[proc_new] = pDict_new
                max_index += 1
        
    # Add processes to param_card
    newFile = os.path.join(outputFolder,'darkcalc_banner.txt')
    with open(paramCard,'r') as f:
        paramCard_data = f.read()
    
    with open(newFile,'w') as f:
        f.write("<DarkCalc version='1.0'>\n")
        f.write("<header>\n")
        f.write("<slha>\n")
        f.write(paramCard_data)
        f.write("</slha>\n")
        f.write("<sigmav>\n")
        # Write header lines
        for pDict in processDict.values():
            header_line = f'#{str(pDict['index']):>{4}},{pDict['name']:<{50}},{pDict['initialPDGs']}_{pDict['finalPDGs']}'
            f.write(header_line+'\n')
        # Now write labels for columns and combine data:
        column_labels = [f"{'x':^10}"]
        data = np.array([])
        for pDict in processDict.values():
            column_labels.append(f"{pDict['index']:^10}")
            if len(data) == 0:
                data = pDict['data']
            else:
                data = np.insert(data,data.shape[-1],
                                 pDict['data'][:,1],axis=1)
        f.write(','.join(column_labels)+'\n')
        # Convert data to string
        data_str = StringIO()
        np.savetxt(data_str, data, fmt='%1.4e', delimiter=',')
        data_str = data_str.getvalue()
        f.write(data_str)
        f.write("</sigmav>\n")
        f.write("</header>\n")
        f.write("</DarkCalc>")

    return newFile


        

