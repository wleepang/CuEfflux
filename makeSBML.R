library(rsbml)

# creates an SBML v2 l1 model from specified model defnition
source('CuEfflux_Func.R')

# set global values
source('CuEfflux_GlobalParams.R')

# after global values are loaded you can change them explicitly here before
# initializing the simulation
modelVersion = 0

# initialize set model parameters and initial state
source('CuEfflux_Init.R')

# define reactions
# produces a list() called 'rxn' that contains a_j and nu_j for each reaction
source('CuEfflux_rxnDef.R')

#nu = set.nu(rxn)
#a  = set.a(rxn)

out = sprintf('CuEffluxModel%s.xml', sub('\\.', 'p', as.character(modelVersion)))

# it appears that sbml can't handle ids with '.' in them
# gsub all '.' in names with '_'
undot = function(n) {gsub('\\.', '_', n)}
redot = function(n) {gsub('_', '\\.', n)}

names(parms) = undot(names(parms))
names(x0) = undot(names(x0))
for (n in names(rxn)) {
  names(rxn[[n]]$nu) = undot(names(rxn[[n]]$nu))
  rxn[[n]]$a = undot(rxn[[n]]$a)
  rxn[[n]]$mods = undot(rxn[[n]]$mods)
}


# create a generic model
sbml = new('SBML')
name(model(sbml)) = sprintf('CuEfflux Model v%.1f', modelVersion)

sbml@model@compartments[['default']] = new('Compartment',
                                           id='default', name='default',
                                           spatialDimensions=3L, size=1, units='litre')
rsbml_write(sbml, out)


# create parameters
for (n in names(parms)) {
  .id = n
  sbml@model@parameters[[.id]] = new('Parameter',
                                   id=.id,
                                   name=n,
                                   value=unname(parms[n]),
                                   constant=T)
}
sbml@model@parameters[['mu']] = new('Parameter', id='mu', name='mu', value=mu, constant=T)
rsbml_write(sbml, out)

# create species
for (n in names(x0)) {
  .id = n
  sbml@model@species[[.id]] = new('Species', 
                                id=.id, 
                                name=n, 
                                compartment='default',
                                initialAmount=unname(x0[n]), 
                                substanceUnits='item', 
                                hasOnlySubstanceUnits=T)
}
rsbml_write(sbml, out)


# create reactions

for (r in 1:length(rxn)) {
  .id = sprintf('r%d', r)
  .name = names(rxn)[r]
  
  cat('rxn id:', .id, 'name:', .name)
  
  if (.name != 'No Reaction Has This Name') {
    speciesList = rxn[[.name]]$nu[which(rxn[[.name]]$nu < 0)]
    reactList = sapply(names(speciesList), function(n){return(new('SpeciesReference', stoichiometry=abs(unname(speciesList[n])), species=n))})
    
    speciesList  = rxn[[.name]]$nu[which(rxn[[.name]]$nu > 0)]
    prodList = sapply(names(speciesList), function(n){return(new('SpeciesReference', stoichiometry=abs(unname(speciesList[n])), species=n))})
    
    speciesList = rxn[[r]]$mods
    modsList = sapply(speciesList, function(n){return(new('ModifierSpeciesReference', id=n, species=n))})
    
    cat(' math:', rxn[[.name]]$a)
    rateExpr  = parse(text=rxn[[.name]]$a, srcfile=NULL)
    #rateParams
    #this gets params based on regex match
    rateParams = names(do.call('c', sapply(c(names(parms), 'mu'), function(p) {grep(p, rxn[[r]]$a, fixed=T)})))
    rateParams = sbml@model@parameters[rateParams]
    
    sbml@model@reactions[[.id]] = new('Reaction', id=.id, name=.name,
                                      reactants=reactList,
                                      products =prodList,
                                      modifiers=modsList,
                                      kineticLaw=new('KineticLaw', math=rateExpr, parameters=rateParams, timeUnits='second', substanceUnits='item'),
                                      reversible=F,
                                      fast=F)
    rsbml_write(sbml, out)
  } else {
    cat()
    cat(' math:', rxn[[.name]]$a, '#### skipped ####')
  }
  cat('\n')
}
#rsbml_write(sbml, out)













