#!/usr/bin/python
"""
# Ivan Domenzain.   Last edited: 2018-10-07
"""
#INPUTS:
ECM_path = '/Users/ivand/Documents/GitHub/ECM_Yeast'
R    = 0.008314472 #Universal constant of gases
Temp = 298.15      #25 celsius
#Load libraries
import os
from os import path
from math import exp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Reaction, ReactionMatcher

#Main Script
eq_api = ComponentContribution(pH=7.5, ionic_strength=0.1)
#Load model data file
os.chdir(ECM_path+'/Models')
data  = pd.read_table('KEGGmodelScaffold.txt')
os.chdir(ECM_path)
i=0
output     = '!!SBtab TableName=Parameter TableType=Quantity Document=S. cervisiae central carbon metabolism SBtabVersion=1.0'+'\n'
output     = output+ '!QuantityType'+'\t'+'!Reaction'+'\t'+'!Compound'+'\t'+'!Value'+'\t'+'!Unit'+'\t'+'!Reaction:Identifiers:kegg.reaction'+'\t'+'!Compound:Identifiers:kegg.compound'+'\t'+'!ID'+'\n'
unit       = 'dimensionless'
compound   = ' '
compoundID = ' '
type       = 'equilibrium constant'

for formula in data.ReactionFormula:
    rxn         = Reaction.parse_formula(formula)
    name        = str(data.ID[i])
    KEGGID      = str(data.KEGGIDs[i])
    reactionStr = name#+'_'+KEGGID    
    parameterID = 'kEQ_'+name
    print(name)
    #Check if reaction is atomically balanced
    if not rxn.check_full_reaction_balancing():
        print('%s is not balanced' % name)
        #raw_input("Press the <ENTER> key to continue...")
        print(' ')
    i = i+1
    #Get standard Gibbs free energy for the rxn
    dG0_prime, dG0_uncertainty = eq_api.dG0_prime(rxn)
    if name == 'RIP1':
        dG0_prime = -38.46 #Edda Klipp's work
    if name == 'RIP1':
        dG0_prime = -38.46 #Edda Klipp's work
    KEQ = exp(-dG0_prime/(R*Temp))
    
    print("dG0' = %.1f uncertainty: %.1f kJ/mol" % (dG0_prime, dG0_uncertainty))
    print ("Keq  = %1f" % (KEQ))   
    #Get reversibility index
    ln_RI  = rxn.reversibility_index()
    print('ln(Reversibility Index) = %.1f\n' % ln_RI)
    output = output+type+'\t'+reactionStr+'\t'+compound+'\t'+str(KEQ)+'\t'+unit+'\t'+ KEGGID+'\t'+compoundID+'\t'+parameterID+'\n'
#Write output:
os.chdir(ECM_path+'/dataFiles')
fid  = open('kEqTable.txt','w')
fid.write(output)
fid.close()
os.chdir(ECM_path)
################################################################################
