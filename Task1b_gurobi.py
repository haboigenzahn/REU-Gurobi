# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:07:01 2017

@author: boigenzahn
"""

from gurobipy import *
import pandas as pd
import numpy as np
import time

starttime = time.time()

file = pd.ExcelFile('Task1b_gurobi.xlsx')

mod = Model('Task1b')

#Compounds
cmpds = ['PW', 'CS', 'EtOH']
feeds = ['PW', 'CS']
prods = ['EtOH']

#Harvest Sites
hs_df = file.parse(skipfooter=677,parse_cols='M')
u_hs = hs_df.values.flatten()
hs = [str(k) for k in u_hs]

#Biorefineries
brfs = ['B1','B2','B3','B4']

#Technologies
techs = ['SHCF','SSCF']

#Time Periods
times = range(1,5)

#----- Parameters ------------
#Availability
#In the form of {(site,compound,time):amount}
alpha_df = file.parse(header=1, skipfooter=677,parse_cols='B:J')
u_alpha = alpha_df.set_index('Sites').T.to_dict('list')
alpha_temp = {str(k): v for k, v in u_alpha.items()}

# Creates tuple keys and assigns values
alpha = {}
for k,v in alpha_temp.iteritems():
    i = 0
    for c in feeds:
        for t in times:
            alpha[(k,c,t)] = v[i]
            i = i+1          
#print alpha

# Demand (liters)
beta = 47250000

# Material loss factor (per period)
gamma = {'PW':0.027,'CS':0.015}

# Annualized capital cost (for a biorefinery producing 2000 tons per day) ($ per year)
zeta = {'SHCF':52780000, 'SSCF':50130000}

# Yield of feedstock (L per ton)
eta = {('PW','SHCF'):283.906,('PW','SSCF'):283.906,('CS','SHCF'):299.048,('CS','SSCF'):299.048}

# Unit inventory cost ($ per ton)
iota = {'PW':1.53,'CS':2.07}

# Fixed unit transportation cost ($ per ton)
kappaf = {'PW':4.29,'CS':10.7}

# Variable unit transportation cost ($ per ton)
kappav = {'PW':0.0737,'CS':0.3059}

# Raw material cost ($ per ton)
llambda = {'PW':65,'CS':59}

# Production cost ($ per ton)
mu = {'SHCF':70,'SSCF':67.7}

# Minimum daily operating capacity (tons)
psil = {'SHCF':1500,'SSCF':1500 }

# Maximum daily operating capacity (tons)
psiu = {'SHCF':3000, 'SSCF':3000}

#Distance to biorefinery 1 in miles
# In the form of (harvest_site, biorefinery):distance
tau_df = file.parse(header=2,skip_footer=677,parse_cols='M:Q')
u_tau = tau_df.set_index('Sites').T.to_dict('list')
temp_tau = {str(k): v for k, v in u_tau.items()}

# Creates (hs,brfs) tuple keys and assigns values
tau = {}
for k,v in temp_tau.iteritems():
    i = 0
    for b in brfs:
        tau[(k,b)] = v[i]
        i=i+1
# print tau

# Working days per period
rho = 85

# ---------- Variables ----------------
# Total annual cost (10^3 $)
TAC = mod.addVar(name="TAC")

# Total cost of feedstock during period t (10^3 $)
Cfeed = mod.addVars(times,name="Cfeed")

# Total cost of transportation during period t (10^3 $)
Ctrans = mod.addVars(times,name="Ctrans")

# Total cost of inventory during period t (10^3 $)
Cinv = mod.addVars(times,name="Cinv")

# Total production cost during period t (10^3 $)
Cprod = mod.addVars(times,name="Cprod")

# Total annual capital cost (10^3 $ per year)
Ccapex = mod.addVar(name="Ccapex")

# Amount of compound sent from harvest site to refinery during period t (10^3 tons)
F = mod.addVars(feeds,hs,brfs,times,name="F")

# Amount of compound harvested during period t (10^3 tons)
H = mod.addVars(feeds,hs,times,name="H")

# Consumption of compound at biorefinery during period (10^3 tons)
G = mod.addVars(feeds,brfs,techs,times,name="G")

# Production of compound at biorefinery during period (10^3 tons)
P = mod.addVars(prods,brfs,times,name="P")

# Daily operating capacity (tons per day)
Q = mod.addVars(brfs,techs,name="Q")

# Inventory of biorefinery at the end of period t (10^3 tons)
Sbr = mod.addVars(feeds,brfs,times,name="Sbr")

# Inventory of harvest site at the end of period t (10^3 tons)
Shs = mod.addVars(feeds,hs,times,name="Shs") 

# Binary selection - will equal 1 if technology m is selected at refinery l
U = mod.addVars(brfs,techs,vtype=GRB.BINARY,name="U")

mod.update()

mod.setObjective(TAC,GRB.MINIMIZE)

# Check - not sure if double sums are working correctly
# Feedstock cost (10^3 $)
for t in times:
    mod.addConstrs(quicksum(H[i,j,t]*llambda[i] for i in feeds for j in hs) == Cfeed[t] for t in times)
# Inventory cost (10^3 $)
for t in times:
    mod.addConstr((quicksum(Sbr[i,l,t]*iota[i] for i in feeds for l in brfs)+quicksum(Shs[i,j,t]*iota[i] for i in feeds for j in hs)) == Cinv[t])
# Inventory cost (10^3 $)    
for t in times:
    mod.addConstr(quicksum(G[i,l,m,t]*mu[m] for i in feeds for l in brfs for m in techs) == Cprod[t])
# Production cost (10^3 $)    
for t in times:
    mod.addConstr(quicksum((kappaf[i]+kappav[i]*tau[(j,l)])*F[i,j,l,t] for i in feeds for j in hs for l in brfs) == Ctrans[t])
# Annualized capital cost based on daily plant capacity
mod.addConstr(quicksum((zeta[m]*Q[l,m]/2000) for m in techs for l in brfs)/1000 == Ccapex)

#Total cost formulation (10^3 $)
mod.addConstr(quicksum(Cfeed[t]+Cinv[t]+Cprod[t]+Ctrans[t] for t in times)+Ccapex == TAC)

# Mass Balance Constraints
# Refinery inventories (10^3 tons)
for i in feeds:
    for l in brfs:
        for t in times:
            if t > 1:
                mod.addConstr((Sbr[i,l,t-1]*(1-gamma[i])+quicksum(F[i,j,l,t] for j in hs)-quicksum(G[i,l,m,t] for m in techs)) == Sbr[i,l,t])
            elif t==1:
                mod.addConstr((Sbr[i,l,t+3]*(1-gamma[i])+quicksum(F[i,j,l,t] for j in hs)-quicksum(G[i,l,m,t] for m in techs)) == Sbr[i,l,t])

# Harvest site inventories (10^3 tons)
for i in feeds:
    for j in hs:
        for t in times:
            if t > 1:
                mod.addConstr((Shs[i,j,t-1]*(1-gamma[i])+H[i,j,t]-quicksum(F[i,j,l,t] for l in brfs)) == Shs[i,j,t])
            elif t==1:
                mod.addConstr((Shs[i,j,t+3]*(1-gamma[i])+H[i,j,t]-quicksum(F[i,j,l,t] for l in brfs)) == Shs[i,j,t])
                
# Cannot harvest more biomass than the site produces (10^3 tons)
mod.addConstrs(alpha[j,i,t] >= H[i,j,t] for i in feeds for j in hs for t in times)      

# Production must meet or exceed demand (liters)
for i in prods:
    for t in times:
        mod.addConstr(quicksum(P[i,l,t] for l in brfs) >= beta)
        
# Amount produced equals amount of feed consumed times yield
for ip in prods:
    for l in brfs:
        for t in times:
            mod.addConstr(quicksum(eta[i,m]*G[i,l,m,t]*1000 for i in feeds for m in techs) == P[ip,l,t])
            
# Capacity constraints
# Amount consumed must be less than the capacity limit of the technology at the biorefinery (tons per day)
mod.addConstrs(quicksum(G[i,l,m,t] for i in feeds)*(1000/rho) <= Q[l,m] for l in brfs for m in techs for t in times)

# Lower capacity of technology (tons)
mod.addConstrs(psil[m]*U[l,m] <= Q[l,m] for l in brfs for m in techs)

# Upper capacity of technology (tons)
mod.addConstrs(psiu[m]*U[l,m] >= Q[l,m] for l in brfs for m in techs)    
       
mod.optimize()

# Print the objective variable
print('\nCost: %g' % mod.objVal)
print("Program Time: %s" % (time.time() - starttime))
# Print the values of variable F
#print(mod.getAttr('x',F))
    
# Print the names and values of all variables
#for v in mod.getVars():
#   print(v.varName, v.x)

# Write python dictionary results to output file
#    target = open('results_v1.txt','w')
#    target.write('\nCost: %g' % mod.objVal)
#    target.write("\n")
#    target.write('Display Variables: %s' % mod.getAttr('x',F))
#    target.close()

# Print results to Excel spreadsheet - further manipulation within the sheet is currently needed to strip punctuation and set delimiters
# F_df = pd.DataFrame.from_dict(mod.getAttr('x',F), 'index')
# writer = pd.ExcelWriter('results_v2_excel.xlsx')
# F_df.to_excel(writer,'Sheet1')
# writer.save()