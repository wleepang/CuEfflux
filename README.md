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

It is assumed that users have a basic understanding of the statistical computing
environment R (freely available for [download](http://www.r-project.org)).

In order to perform simulations users must first install the `deSolve` package:
```{r}
install.packages('deSolve')
```

## Basic Simulation
For an example of basic simulation code, refer to the file:
 - odeCuEfflux.R

Before beginning, utility functions must be imported from:
 - CuEfflux_Func.R

which provides:
 - `dxdt()` : the differential step function
 - `.T()` : for the calculation of 'total' abunance of a species
 - `setStoic()` : sets individual reaction stoichiometries
 - `set.a()`, `set.nu()` : creates the propensity/rate vector `a`, and stoichiometry
   matrix `nu` from reaction definitions (see below)

The model is defined as reaction network that is simulated as a system of
ordinary differential equations.

This definition is separated into three files:
 - CuEfflux_GlobalParams.R
 - CuEfflux_Init.R
 - CuEfflux_rxnDef.R

which are imported into the simulation code using `source()`.

**GlobalParams**: specifies parameters used throughout the model.  This includes all
default reaction specific rate constants, and boolean flags that determine the
type of genetic background being simulated.  After import, a named vector called
`params` will be added to the workspace.

**Init**: specifies initial conditions for all species.  This file also initializes
parameter values based on the 'modelName' parameter defined in GlobalParams.
For example, if a reaction is not used in a model, its reaction rate is set to a
value of 0.  After import a named vector called `x0` will be added to the
workspace and the vector `params` will have been modified accordingly.

**rxnDef**: specifies all reactions used by all models.  After import, a named list
called `rxn` will be added to the workspace.  Each element of `rxn` will have
the following structure:

```{r}
rxn$ReactionName
 .. $ a  : chr()
 .. $ nu : num()
```

where
 - `$a` is a character string defining the reaction rate law
 - `$nu` is a numeric vector whose names are the species involved in the reaction
   and values are their respective stoichiometries.  The stoichiometry
   convention is reactants < 0, products > 0.

After reactions are defined, the rate vector `a` and stoichiometry matrix `nu`
are created by calling:
```{r}
a = set.a(rxn)
nu = set.nu(rxn)
```

The elements of the rate vector `a` are the defined rate laws from `rxn` and
retain the reaction names defined in `rxn`.

The rows and columns of `nu` correspond to species and reactions, respectively,
and are named accordingly.

Simulations are carried out via a call to `ode()` (from the deSolve package):
```{r}
out = ode(x0, 
          seq(0,18000,by=100), 
          dxdt, 
          c(list(nu=nu, a=a), parms), 
          method='daspk')
```

Here, the simulation is performed using initial conditions `x0` from timepoints
0 to 18000 seconds, advancing every 100 seconds.  The differential step function
`dxdt()` is called with additional parameters `c(list(nu=nu, a=a), params)`.  The
simulation algorithm is set to `daspk` which is best suited for stiff ODE
systems.

The result `out` is a matrix-like object with rows and columns as timepoints and
species, respectively.

To plot the simulated time profile of the species `Cu`, the intracellular copper
abundance, first a time vector in units of minutes is generated:
```{r}
t.abs = out[,'time'] / 60
```

Then the raw `Cu` trajectory is plotted via:
```{r}
plot(t.abs, out[,'Cu'], type='l')
```

Here it is important to note that `Cu` exists as both an isolated species and
bound to other molecules such as metallochaperones and export pumps.  Thus, it
is necessary to 'totalize' all `Cu` species to get the correct intracellular
abundance:
```{r}
plot(t.abs, .T('Cu', out), type='l')
```

Finally, the model accounts for cellular growth.  To get the intracellular `Cu`
abundance on a per cell basis (assuming a volume of 1 is one cell):
```{r}
plot(t.abs, .T('Cu', out)/out[, 'V'], type='l')
```

## Summary
The sequence of steps above describe the core simulation workflow.  All other
analyses - e.g. parameter sweeps, strain comparisons, sensitivity analysis - 
merely extend this workflow.

## Appendix 1: SBML Export
Code is available to convert the model as it exists in R to [SBML](www.sbml.org)
a commonly used format for easy exchange and simulation of models in a variety
of computational tools.

### Requirements
 - Cu Efflux Model R Code
 - R 2.11 (or higher)
 - R Packages:
   - rsbml

### Coversion
To convert the R model to an SBML file simply run:
```{r}
source('makeSBML.R')
```

This will produce a file called `CuEfflux.xml` that is SBML L2V1 compliant which
can then be imported and further modified/tested using SBML aware software tools
such as:
 - COPASI
 - CellDesigner

#### Caveats
If any modifications to the model are made one must take the following into
consideration prior to conversion to SBML.

 1. rsbml does not check SBML consistency. An SBML file may be produced but may
    generate errors during import.
 2. Mathematical expressions for reaction KineticLaws converted from R to
    SBML/MathML must be in fully expanded infix form - e.g.:
```
a*(b + c)
```
    will produce libSBML errors, whereas
```
a*b + a*c
```
    will covert to SBML/MathML.  This limitation is not specified anywhere in
    either the rsbml, SBML, or libSBML documentation.
 3. Object names cannot have '.'.  Substituting with '_' passes.