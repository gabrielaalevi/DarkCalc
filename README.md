DarkCalc is a Linux-based numerical tool to solve the Boltzmann equations for a generic model. We use MadGraph and MadDM to compute the thermally-averaged cross sections and decay widths relevant for the model.

## Installation

Prior to installation, the user must pre-install the packages `NumPy`, `SciPy`, `six`, `pyslha` and `git`, all compatible with Python 3. Firstly, it is necessary to clone the GitHub repository, via the command line:
```
mydir$ git clone https://github.com/gabrielaalevi/DarkCalc.git
```
then navigate to the DarkCalc directory and run the installer script:
```
mydir$ cd DarkCalc
DarkCalc$ ./installer.sh
```
This will prompt a message:
```
Install MadGraph (y/n)?
```
By pressing `y` and `Enter`, the installation will be completed.

## Input

There are two essential inputs for DarkCalc: the UFO model files and the input file. 

The UFO model must be pre-downloaded by the user, from the FeynRules Model database. DarkCalc comes with a `models` directory to store such files, but the user can personalized the path to the UFO model in the input file, allowing for storage of models in arbitrary locations. By default, DarkCalc comes with the `DMSimp-t` model pre-installed.

The physical parameters and execution settings are specified via the `input_parameters.ini` file, within the root directory `DarkCalc`. A default example file is distributed with DarkCalc, to serve as a template. The available options in the input file are:

*Options*: definition of paths, and options for turning on and off some DarkCalc features.

- `outputFolder` (path): defines the directory in which DarkCalc will store all output files related to the current run, excluding the input file itself and the solution file.
  
- `MadGraphPath` (path): specifies the installation path for MadGraph. If the user executed the installer script, it should be set to \texttt{MG5}.
  
- `skipMadDM` (True/False): if False, DarkCalc will execute MadDM to compute thermally-averaged cross sections and decay widths. If True, the code does not execute MadDM, and instead searches for cross section and decay width values inside the `outputFolder`. This option is useful for re-running only the Boltzmann equation integration, provided that particle masses and couplings remain unchanged.
  
- `addConversion` (True/False): controls the inclusion of conversion processes in the Boltzmann equation. if True, such interactions are included in the solution. If False, these reactions are omitted.
  
*Model*: definition of the model characteristics, such as the UFO model and the particles present in the simulation.

- `modelDir` (path): path to the UFO model directory to be used in the simulation.
  
- `darkmatter` (string): label of the Dark Matter particle.
  
- `bsmParticles` (list of strings): list of all Beyond-the-Standard-Model particles present in the model, identified by labels.
  
- `computeWidths` (list of strings): list of BSM particles whose decay width should be computed.
  
*SetParameters*: specification of the model parameters, such as masses and coupling values.

- mass parameters: particle masses are set via lines of the form

m(label) = (value in GeV)

where (label) identifies the particle.

- other parameters: any additional model-specific parameters (such as coupling strengths), are defined by:
  
(parameter name) = (value)

  A complete list of available parameters can be found in the UFO model documentation on the FeynRules database.
  
*SolverParameters*: options for the numerical integrator and for output formatting.

- `initialConditions` (dictionary): a Python dictionary specifying the initial conditions for each particle's yield, such as:
  
  `initialConditions` = {(label): *condition*, ...}
  
where 'condition' can be *thermal* (for thermal equilibrium) or a numeric value setting the desired yield.

- `atol`, `rtol` (floats): absolute and relative tolerances for the numerical integration of the Boltzmann equations. Defaults are both 10^{-10}, if unspecified.
  
- `T0`, `Tf` (floats): initial and final temperatures for integration, in GeV.
  
- `method` (string): specifies the numerical integration method. Acceptable values are those supported by `SciPy`s function `solve_ivp`. Default is 'Radau'.
  
- `nsteps` (integer): number of evaluation steps for the numerical solver.
  
- `outputFile` (string): name for the output file, along with the file type identifier suffix. By default, the file is saved in the root directory `DarkCalc`.
  
- `extendedOutput` (True/False): if True, DarkCalc saves the contribution from each collision term, for each interaction, in addition to the solution of the Boltzmann equations.
  
## Execution and Output

To execute the code, the user must run in the command line:

```DarkCalc$ python3 ./main.py```

This command will prompt the execution of DarkCalc. The output file will be saved by name `outputFile`, in the format chosen by the user in this variable, inside the root `DarkCalc` directory. To plot the solution, DarkCalc comes with a Jupyter notebook named `plot.ipynb`, inside the `examples` folder.
