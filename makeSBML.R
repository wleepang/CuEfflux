library(rsbml)

# creates an SBML v2 l1 model from specified model defnition
source('CuEfflux_Func.R')

# set global values
source('CuEfflux_GlobalParams.R')

# after global values are loaded you can change them explicitly here before
# initializing the simulation
modelVersion = 1.3

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

nu = set.nu(rxn)
a  = set.a(rxn)

# create a generic model
sbml = new('SBML')

# create parameters
for (n in names(parms)) {
  sbml@model@parameters[[n]] = new('Parameter',
                                   id=n,
                                   name=n,
                                   value=unname(parms[n]),
                                   constant=T)
}
sbml@model@parameters[['mu']] = new('Parameter', id='mu', name='mu', value=mu, constant=T)

# create species
for (n in names(x0)) {
  sbml@model@species[[n]] = new('Species', 
                                id=n, 
                                name=n, 
                                initialAmount=unname(x0[n]), 
                                substanceUnits='molecules', 
                                hasOnlySubstanceUnits=T)
}

# create reactions

reactList = rxn[[name]]$nu[which(rxn[[name]]$nu < 0)]
prodList  = rxn[[name]]$nu[which(rxn[[name]]$nu > 0)]
rateExpr  = rxn[[name]]$a
#rateParams
#this gets params based on regex match
do.call('c', sapply(names(parms), function(p) {grep(p, rxn[[1]]$a, fixed=T)}))


new('Reaction', id=id, name=name,
    reactants=reactList,
    products =prodList,
    kineticLaw=new('KineticLaw', math=rateExpr, parameters=rateParams),
    reversible=F
    fast=F)














