import os
from tools.logger import logger
import tempfile
import subprocess
import numpy as np



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

    if 'addConversion' in parser['Options']:
        addConversion = bool(parser['Options']['addConversion'])
    else:
        addConversion = False

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
        return mergeOutput(outputFolder, addConversion)

def mergeOutput(outputFolder : str, addConversion: bool) -> str:
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

    lines = [l.strip() for l in sigmaV.splitlines() if l.strip()]
    # Get process dictionary        
    comment_lines =  [l for l in lines if l[0] == '#']
    proc_lines = [l[1:].strip() for l in comment_lines if (l.count('#') == 1)]
    processDict = {}
    for l in proc_lines:    
        proc_index,proc_name,proc_pdgs = (l.split(',',2))
        proc_pdgs = proc_pdgs.replace('[','').replace(']','').replace('-', '')
        iPDGs,fPDGs = proc_pdgs.split('_')
        initialPDGs = list(map(int, iPDGs.split(',')))
        finalPDGs = list(map(int, fPDGs.split(',')))
        count_reaction = 0
        for j in processDict.items():
            dict_initial_pdgs = j[1]['initialPDGs']
            dict_final_pdgs = j[1]['finalPDGs']
            if (int(initialPDGs[0]) == int(dict_initial_pdgs[0])) and float(initialPDGs[1]) == float(dict_initial_pdgs[1]) and float(finalPDGs[0]) < 30 and float(finalPDGs[1]) < 30 and float(dict_final_pdgs[0])<30 and float(dict_final_pdgs[1])<30:
                count_reaction += 1
        if count_reaction == 0:
            processDict[proc_index] = {'name' : proc_name.strip(), 
                                    'initialPDGs' : initialPDGs, 
                                    'finalPDGs' : finalPDGs}
            processes_data = np.genfromtxt(lines, delimiter=',', skip_header=len(comment_lines), 
                                    names=True)
        for proc_index,pInfo in processDict.items():
            if addConversion == True:
                counter = 0 #we only need one co-annihilation interaction to calculate the conversion taacs, so we set a counter
                if abs(pInfo['initialPDGs'][0]) != abs(pInfo['initialPDGs'][1]) and counter == 0:
                    counter = counter + 1
                    conv_index = len(proc_lines) + 1
                    conv_line = '#   ' + str(conv_index) + ',conversion                                        ,['+str(pInfo['initialPDGs'][0])+', '+str(pInfo['finalPDGs'][0])+']_['+str(pInfo['initialPDGs'][1])+', '+str(pInfo['finalPDGs'][1])+']\n'
                    conv_taacs = processes_data[proc_index]
                    processes =  ''
                    taacs_header = 'x,'
                    for k in processDict.items():
                        processes = processes + '#   ' + str(k[0]) + ',' + k[1]['name'] + '                                          ,[' + str(k[1]['initialPDGs'][0])+',' + str(k[1]['initialPDGs'][1]) +']_[' + str(k[1]['finalPDGs'][0]) + ',' + str(k[1]['finalPDGs'][1]) + ']\n'
                        taacs_header = taacs_header + ' ' + str(k[0]) + ','
                    processes = processes + conv_line
                    taacs_header = taacs_header + ' ' + str(conv_index) + '\n'
                    original_processes, taacs = sigmaV.split('x,  1', 1)
                    header_original, taacs_values = taacs.split(str(len(proc_lines)), 1)
                    taacs_lines = [l.strip() for l in taacs_values.splitlines() if l.strip()]
                    new_taacs_lines = ''
                    for i in range(0, len(taacs_lines)):
                        new_taacs_lines = new_taacs_lines + str(processes_data[i][0])
                        for m in processDict.items():
                            new_taacs_lines = new_taacs_lines + ',' + str(processes_data[i][m[0]])
                        new_taacs_lines = new_taacs_lines + ',' + str(conv_taacs[i]) + '\n'
            if addConversion == False:
                processes =  ''
                taacs_header = 'x'
                for k in processDict.items():
                    processes = processes + '#   ' + str(k[0]) + ',' + k[1]['name'] + '                                          ,[' + str(k[1]['initialPDGs'][0])+',' + str(k[1]['initialPDGs'][1]) +']_[' + str(k[1]['finalPDGs'][0]) + ',' + str(k[1]['finalPDGs'][1]) + ']\n'
                    taacs_header = taacs_header + ', ' + str(k[0])
                taacs_header = taacs_header + '\n'
                original_processes, taacs = sigmaV.split('x,  1', 1)
                header_original, taacs_values = taacs.split(str(len(proc_lines)), 1)
                taacs_lines = [l.strip() for l in taacs_values.splitlines() if l.strip()]
                new_taacs_lines = ''
                for i in range(0, len(taacs_lines)):
                    new_taacs_lines = new_taacs_lines + str(processes_data[i][0])
                    for m in processDict.items():
                        new_taacs_lines = new_taacs_lines + ',' + str(processes_data[i][m[0]])
                    new_taacs_lines = new_taacs_lines + '\n'

    newFile = os.path.join(outputFolder,'darkcalc_banner.txt')
    with open(newFile,'w') as f:
        f.write("<DarkCalc version='1.0'>\n")
        f.write("<header>\n")
        f.write("<slha>\n")
        f.write(paramCard)
        f.write("</slha>\n")
        f.write("<sigmav>\n")
        f.write(processes)
        f.write(taacs_header)
        f.write(new_taacs_lines)
        f.write("</sigmav>\n")
        f.write("</header>\n")
        f.write("</DarkCalc>")

    return newFile


        

