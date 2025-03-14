import os
from tools.logger import logger
import tempfile
import subprocess



def runMadDM(parser : dict) -> str:
    """
    Run MadDM to compute widths and relevant collision cross-sections.
    
    :param parser: Dictionary with parser sections.
    
    :return: True if successful. Otherwise False.
    """
        
    #Get run folder:    
    outputFolder = os.path.abspath(parser['Options']['outputFolder'])
    if not os.path.isdir(outputFolder):
        os.makedirs(outputFolder)    
    
    modelDir = os.path.abspath(parser['Model']['modelDir'])
    if not os.path.isdir(modelDir):
        # Check for model restriction in the name
        mDir = modelDir.rsplit('-',1)[0]
        if not os.path.isdir(mDir):
            logger.error(f'Model folder {modelDir} (or {mDir}) not found')
            return False

    dm = parser['Model']['darkmatter']
    if 'bsmParticles' in parser['Model']:
        bsmList = str(parser['Model']['bsmParticles']).split(',')
        # Make sure the dmPDG does not appear twice
        bsmList = [p for p in bsmList[:] if p != dm]
    else:
        bsmList = []

    #Generate commands file:       
    commandsFile,cFilePath = tempfile.mkstemp(suffix='.txt', prefix='maddm_commands_', dir=outputFolder)    
    os.close(commandsFile)
    commandsFileF = open(cFilePath,'w')
    commandsFileF.write(f'import model {modelDir}\n')
    commandsFileF.write(f'define darkmatter {dm}\n')
    for p in bsmList:
        commandsFileF.write(f'define coannihilator {p}\n')
    commandsFileF.write('generate relic_density\n')
    commandsFileF.write(f'output {outputFolder}\n')
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
        return False
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
    sigmaVFile = os.path.join(outputFolder,'output','taacs.csv')
    if not os.path.isfile(sigmaVFile):
        logger.error(f"Error computing sigmaV with MadDM ({sigmaVFile} not found)")
        return None
    else:
        return mergeOutput(outputFolder)

def mergeOutput(outputFolder : str) -> str:
    """
    Combines the param_card and the sigmaVFile into a single file,
    similar to the MadGraph banner.
    Returns the new file name
    """

    paramCard = os.path.join(outputFolder,'Cards','param_card.dat')
    if  not os.path.isfile(paramCard):
        logger.error(f"Param card ({paramCard} not found)")
        return None
    sigmaVFile = os.path.join(outputFolder,'output','taacs.csv')
    if not os.path.isfile(sigmaVFile):
        logger.error(f"SigmaV file ({sigmaVFile} not found)")
        return None

    with open(paramCard,'r') as f:
        paramCard = f.read()
    with open(sigmaVFile,'r') as f:
        sigmaV = f.read()

    newFile = os.path.join(outputFolder,'darkcalc_banner.txt')
    with open(newFile,'w') as f:
        f.write("<DarkCalc version='1.0'>\n")
        f.write("<header>\n")
        f.write("<slha>\n")
        f.write(paramCard)
        f.write("</slha>\n")
        f.write("<sigmav>\n")
        f.write(sigmaV)
        f.write("</sigmav>\n")
        f.write("</header>\n")
        f.write("</DarkCalc>")

    return newFile


        

