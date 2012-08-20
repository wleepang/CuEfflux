# Model code for Cu Efflux in Halobacterium salinarum

## Developed by
Wyming Lee Pang, PhD
at the
Institute for Systems Biology

## Overview
The following is code that defines a quantitative model of copper efflux in the
halophillic archaea _Halobacterium salinarum_, with specific focus on how the
dynamics of gene expression of efflux machinery and intracellular copper levels
change in response to varying levels of metallochaperone expression.

A publication that for which this model was created and provides more scientific
detail is currently under review at PLoS Computational Biology.

### Requirements:
 - R (2.11 or higher)
 - R packages
    - deSolve
    - matlab (optional)

It is assumed that users have a basic understanding of the statistical computing
environment R (freely available for download at [http://www.r-project.org]).

In order to perform simulations users must first install the 'deSolve' package:
> install.packages('deSolve')

Convenient timing functions tic() and toc() are provided by the 'matlab'
package, whose installation is optional.

## Basic Simulation
For an example of basic simulation code, refer to the file:
 - odeCuEfflux.R

Before beginning, utility functions must be imported from:
 - CuEfflux_Func.R

which provides:
 - dxdt() : the differential step function
 - .T() : for the calculation of 'total' abunance of a species
 - setStoic() : sets individual reaction stoichiometries
 - set.a(), set.nu() : creates the propensity/rate vector a, and stoichiometry
   matrix nu from reaction definitions (see below)

The model is defined as reaction network that is simulated as a system of
ordinary differential equations.

This definition is separated into three files:
 - CuEfflux_GlobalParams.R
 - CuEfflux_Init.R
 - CuEfflux_rxnDef.R

which are imported into the simulation code using source().

GlobalParams: specifies parameters used throughout the model.  This includes all
default reaction specific rate constants, and boolean flags that determine the
type of genetic background being simulated.  After import, a named vector called
'params' will be added to the workspace.

Init: specifies initial conditions for all species.  This file also initializes
parameter values based on the 'modelName' parameter defined in GlobalParams.
For example, if a reaction is not used in a model, its reaction rate is set to a
value of 0.  After import a named vector called 'x0' will be added to the
workspace and the vector 'params' will have been modified accordingly.

rxnDef: specifies all reactions used by all models.  After import, a named list
called 'rxn' will be added to the workspace.  Each element of 'rxn' will have
the following structure:

rxn$<name>
 .. $ a  : chr()
 .. $ nu : num()

where
 - $a is a character string defining the reaction rate law
 - $nu is a numeric vector whose names are the species involved in the reaction
   and values are their respective stoichiometries.  The stoichiometry
   convention is reactants < 0, products > 0.

After reactions are defined, the rate vector 'a' and stoichiometry matrix 'nu'
are created by calling:
> a = set.a(rxn)
> nu = set.nu(rxn)

The elements of the rate vector 'a' are the defined rate laws from 'rxn' and
retain the reaction names defined in 'rxn'.

The rows and columns of 'nu' correspond to species and reactions, respectively,
and are named accordingly.

Simulations are carried out via a call to ode() (from the deSolve package):
> out = ode(x0, 
            seq(0,18000,by=100), 
            dxdt, 
            c(list(nu=nu, a=a), parms), 
            method='daspk')

Here, the simulation is performed using initial conditions 'x0' from timepoints
0 to 18000 seconds, advancing every 100 seconds.  The differential step function
dxdt() is called with additional parameters c(list(nu=nu, a=a), params).  The
simulation algorithm is set to 'daspk' which is best suited for stiff ODE
systems.

The result 'out' is a matrix-like object with rows and columns as timepoints and
species, respectively.

To plot the simulated time profile of the species 'Cu', the intracellular copper
abundance, first a time vector in units of minutes is generated:
> t.abs = out[,'time'] / 60

Then the raw 'Cu' trajectory is plotted via:
> plot(t.abs, out[,'Cu'], type='l')

Here it is important to note that 'Cu' exists as both an isolated species and
bound to other molecules such as metallochaperones and export pumps.  Thus, it
is necessary to 'totalize' all 'Cu' species to get the correct intracellular
abundance:

> plot(t.abs, .T('Cu', out), type='l')

Finally, the model accounts for cellular growth.  To get the intracellular 'Cu'
abundance on a per cell basis (assuming a volume of 1 is one cell):

> plot(t.abs, .T('Cu', out)/out[, 'V'], type='l')

## Summary
The sequence of steps above describe the core simulation workflow.  All other
analyses - e.g. parameter sweeps, strain comparisons, sensitivity analysis - 
merely extend this workflow.