#Cory Williams


# IMPORT A MODEL

from __future__ import print_function
import cobra
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import itertools
import os 

model = cobra.io.load_json_model("Recon3D.json")


# Do single gene deletions
single_del = cobra.flux_analysis.single_reaction_deletion(model)
single_gene_del = cobra.flux_analysis.single_gene_deletion(model)

# Useful info about data
print("How many reactions are in the model:", len(model.reactions))
print("How many metabolites are in the model:",len(model.metabolites))
print("How many genes are in the model:",len(model.genes))
print(model.objective)

#reaction in model that contains ACE/ACE2 genes
model.reactions.get_by_id('RE2445E')
#ACE2 gene
model.genes.get_by_id('59272_AT1')
#ACE gene
model.genes.get_by_id('1636_AT1')

#assigning results script to variable
solution = model.optimize()

#knocking out all associated reactions
rint('Reaction Knockouts:')
print('No knock out: ', model.optimize())
print()
with model:
    model.reactions.RE2445E.knock_out()
    print('Reaction: RE2445E knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.GLNt4.knock_out()
    print('Reaction: GLNt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.PROt4.knock_out()
    print('Reaction: PROt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.ALAt4.knock_out()
    print('Reaction: ALAt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.ASNt4.knock_out()
    print('Reaction: ASNt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.THRt4.knock_out()
    print('Reaction: THRt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.RE0938E.knock_out()
    print('Reaction: RE0938E knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.SERt4.knock_out()
    print('Reaction: SERt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.LEUt4.knock_out()
    print('Reaction: LEUt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.PHEt4.knock_out()
    print('Reaction: PHEt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.GLYt4.knock_out()
    print('Reaction: GLYt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.ILEt4.knock_out()
    print('Reaction: ILEt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.TYRt4.knock_out()
    print('Reaction: TYRt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.VALt4.knock_out()
    print('Reaction: VALt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.METt4.knock_out()
    print('Reaction: METt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.RE0937E.knock_out()
    print('Reaction: RE0937E knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.RE0936E.knock_out()
    print('Reaction: RE0936E knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.TRPt4.knock_out()
    print('Reaction: TRPt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.RE2445E.knock_out()
    print('Reaction: RE2445E knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.CYSt4.knock_out()
    print('Reaction: CYSt4 knocked out: ', model.optimize())
    print() 
with model:
    model.reactions.SELMETHte.knock_out()
    print('Reaction: SELMETHte knocked out: ', model.optimize())
    print() 

#creating gene list and knocking out all genes in list
genelist = ['59272_AT1','1636_AT1','81539_AT1','11254_AT1','54407_AT1', '55089_AT1', '340024_AT1', '57393_AT1']
for i in genelist:
    with model:
        print()
        model.genes.get_by_id(i).knock_out()
        print('Gene:',i, 'knocked out: ', model.optimize())
        
#calculating flux expressions  
model.reactions.RE2445E.flux_expression

      
        same_flux = model.problem.Constraint(
    model.reactions.RE2445E.flux_expression - model.reactions.GLNt4.flux_expression,
    lb=-1000,
    ub=0)
model.add_cons_vars(same_flux)


solution = model.optimize()
print(solution.fluxes['RE2445E'], solution.fluxes['GLNt4'],
      solution.objective_value)