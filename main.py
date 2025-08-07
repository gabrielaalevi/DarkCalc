#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys,os
from tools.configParserWrapper import ConfigParserExt
from tools.logger import logger,setLogLevel
from tools.maddm_interface import runMadDM
from tools.modelData import ModelData
from boltz.boltzSolver import runSolver
from boltz.boltzmannEq import computeCollisionTerms,computeDecayTerms
import thermal.equilibriumDensities as eqDensitities
import multiprocessing
import time,datetime
import numpy as np




def saveSolutions(parser : dict, solution, model : ModelData) -> bool:

    pars = parser['SolverParameters']
    extended = bool(pars['extendedOutput'])
    compDict = model.componentsDict
    labels = [comp.label for comp in sorted(compDict.values(), key= lambda comp: comp.ID) if comp.ID != 0]
    x_sol = solution.t
    y_sol = solution.y
    headerList = ['x'] + [f'Y({label})' for label in labels] 
    data = np.array(list(zip(x_sol,*y_sol[1:])))
    Ytot = sum(data[-1][1:])
    for i,label in enumerate(labels):
        Y_i = data[-1][i+1]
        om_i = 0.12*Y_i/(6.8e-13)
        logger.info(f"Omega*h^2({label}) = {om_i:1.3e}")
        
    Ytot = sum(data[-1][1:])
    omh2 = 0.12*Ytot/(6.8e-13)
    logger.info(f'\n\nOmega*h^2 = {omh2:1.4g}\n')

    if extended:
        mDM = compDict[model.dmPDG].mass
        Di_vec = []
        Cij_vec = []
        proc_names = None
        for i,x in enumerate(x_sol):
            H = eqDensitities.H(mDM/x) #hubble rate at temperature T
            dsdx = np.abs(eqDensitities.dSdx(x, mDM)) #variation of entropy with x
            Y = y_sol[:,i]
            # Print decay terms
            Di_vec.append([(dec['decay']/(3*H))*dsdx for dec in computeDecayTerms(x,Y,model)][1:])
            c = computeCollisionTerms(x,Y,model)
            if proc_names is None:                
                proc_names = [name for cc in c for name in cc.keys()]
            proc_values = [(sigma/(3*H))*dsdx 
                           for cc in c for sigma in cc.values()]
            # Print collision terms for each process
            Cij_vec.append(proc_values)
        
        headerList += [f'D({labelA})*(dsdx/3H)' for labelA in labels]
        headerList += [f'C({proc_name})*(dsdx/3H)' for proc_name in proc_names]
        
        data = np.hstack((data,Di_vec))
        data = np.hstack((data,Cij_vec))
        
    header = ','.join(headerList)
    outFile = os.path.abspath(pars['outputFile'])
    np.savetxt(outFile,data,header=header,fmt='%1.4e',delimiter=',')  

    return True

def runSolution(parser : dict) -> bool:
    """
    Run MadDM, load the model and solve the Boltzmann equations

    :param parser: Dictionary with parser sections.
    
    :return: Dictionary with run info. False if failed.
    """

    if 'skipMadDM' in parser['Options']:
        skipMadDM = bool(parser['Options']['skipMadDM'])
    else:
        skipMadDM = False
    if not skipMadDM:
        logger.info("Running MadDM")
        outputFile = runMadDM(parser)
        logger.info("Finished MadDM run")
    else:
        outputFolder = os.path.abspath(parser['Options']['outputFolder'])
        outputFile = os.path.join(outputFolder,'darkcalc_banner.txt')
    logger.info("Loading model")
    model = ModelData.loadModel(parser, outputFile)
    logger.info("Model loaded")
    t0 = time.time()
    logger.info("Solving Boltzmann equations")
    sol = runSolver(parser,model)
    if sol.success:
        dt = (time.time()-t0)
        logger.info(f"Solved Boltzmann equations in {dt:1.1f} s")
        logger.info("Saving solutions")
        saveSolutions(parser,sol,model)
    else:
        logger.error(f"Error solving Boltzmann equations:\n {sol.message}\n")
    return sol


def main(parfile,verbose):
   
    level = verbose
    setLogLevel(level)    

    parser = ConfigParserExt(inline_comment_prefixes="#")   
    # Define default solver parameters
    ret = parser.read_dict({'SolverParameters' : {'atol' : 0.0, 
                                                  'rtol' : 1e-4, 
                                                  'method' : 'Radau', 
                                                  'nsteps' : 100}})
    ret = parser.read(parfile)
    if ret == []:
        logger.error( "No such file or directory: '%s'" % args.parfile)
        sys.exit()
            
    #Get a list of parsers (in case loops have been defined)    
    parserList = parser.expandLoops()

    # Start multiprocessing pool
    ncpus = -1
    if parser.has_option("options","ncpu"):
        ncpus = int(parser.get("options","ncpu"))
    if ncpus  < 0:
        ncpus =  multiprocessing.cpu_count()
    ncpus = min(ncpus,len(parserList))
    pool = multiprocessing.Pool(processes=ncpus)
    if ncpus > 1:
        logger.info('Running %i jobs in parallel with %i processes' %(len(parserList),ncpus))
    else:
        logger.info('Running %i jobs in series with a single process' %(len(parserList)))

    now = datetime.datetime.now()
    children = []
    for irun,newParser in enumerate(parserList):
        # Create temporary folder names if running in parallel
        parserDict = newParser.toDict(raw=False)
        logger.debug('submitting with pars:\n %s \n' %parserDict)
        p = pool.apply_async(runSolution, args=(parserDict,),)
        children.append(p)

#     Wait for jobs to finish:
    output = [p.get() for p in children]
    logger.info("Finished all runs (%i) at %s" %(len(parserList),now.strftime("%Y-%m-%d %H:%M")))

    return output


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Runs MadDM to compute cross-sections and solve the Boltzmann equations for the parameters defined in the parameters file." )
    ap.add_argument('-p', '--parfile', default='input_parameters.ini',
            help='path to the parameters file [input_parameters.ini].')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')



    t0 = time.time()

    args = ap.parse_args()
    output = main(args.parfile,args.verbose)
            
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
