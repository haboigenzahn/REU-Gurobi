# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:28:15 2017

@author: boigenzahn
"""
from gurobipy import *
import pandas as pd
import numpy as np
from collections import defaultdict

file = pd.ExcelFile('Task2a_gurobi.xlsx')

mod = Model('Task2a')

#Compounds
cmpds = ['CS', 'pellet','EtOH']
feeds = ['CS']
itrmds = ['pellet']
prods = ['EtOH']

#Harvest Sites
hs_df = file.parse(sheetname='Sheet1', header=None, skipfooter = 1, parse_cols='B')
u_hs = hs_df.values.flatten()
hs = [str(k) for k in u_hs]
#print hs

#Depots
dpt_df = file.parse(sheetname='Sheet1',header=0,parse_cols='U')
u_dpt = dpt_df.values.flatten()
dpts = [str(k) for k in u_dpt]
#print dpts

#Biorefineries
brfs = ['B29']

#Technologies
techs = ['pre','SSCF']
brfs_techs = ['SSCF']
dpts_techs=['pre']

#Time Periods
times = range(1,5)

#----- Parameters ------------
#Availability
#In the form of {(site,compound,time):amount}
alpha_df = file.parse(header = 1, sheetname='Sheet1',parse_cols='B:F')
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

# Demand (10^6 liters)
beta = 210

# Material loss factor (per period)
gamma = {'CS':0.015, 'pellet':0.006, 'EtOH':0}

# Annualized capital cost (for a biorefinery producing 2000 tons per day) (10^6$ per year)
zeta = {'pre':0.8923, 'SSCF':50.13}

# Yield of compound i to compound ii using technology m (ton/ton or L/ton)
# Conversion to defaultdict makes all key values not found in the dictionary equal to 0
eta = {('CS', 'EtOH','SSCF'):299.048,('CS','pellet','pre'):0.99,('pellet','EtOH','SSCF'):299.048}
eta = defaultdict(lambda: 0, eta)

# Unit inventory cost ($ per ton)
iota = {'CS':1.433,'pellet':0.667}

# Fixed unit trucking cost ($ per ton)
kappaTf = {'CS':6.77}

# Variable unit trucking cost ($ per ton)
kappaTv = {'CS':0.319}

# Fixed unit rail cost ($ per ton)
kappaRf = {'pellet':17.34}

# Variable unit rail cost ($ per ton)
kappaRv = {'pellet':0.018}

# Raw material cost ($ per ton)
llambda = {'CS':40}

# Production cost ($ per ton)
mu = {'pre':15.15,'SSCF':67.7}

# Minimum daily operating capacity (tons)
psil = {'pre':100,'SSCF':800 }

# Maximum daily operating capacity (tons)
psiu = {'pre':3000, 'SSCF':9000}

#TODO: Check that the fact that these columns have the same name isn't causing the column to actually be duplicated/excluded
tau_df = file.parse(sheetName='Sheet1',skip_footer=1,parse_cols='Q,R,U,V')
#Distance from harvest site to biorefinery (miles)
u_taujl = tau_df.iloc[: , 0:2].set_index('Sites').T.to_dict('list')
taujl = {str(k): v[0] for k, v in u_taujl.items()}
#print taujl
#Distance from depot to biorefinery (miles)
u_taukl = tau_df.iloc[: , 2:4].set_index('Depots').T.to_dict('list')
taukl = {str(k): v[0] for k, v in u_taukl.items()}
#print taukl

tau_df2 = file.parse(sheetname=1,header=2,parse_cols='C:EH')
u_taujk = tau_df2.T.to_dict('list')
temp_taujk = {str(k):v for k,v in u_taujk.items()}
# Creates (hs,dpts) tuple keys and assigns values
taujk = {}
for h,v in temp_taujk.iteritems():
    i = 0
    for k in dpts:
        taujk[(h,k)] = v[i]
        i=i+1
#print taujk        

# Working days per period
rho = 85

# --------------------- More Sets -------------------
# TODO: check that these re-assignments are working and not butchering the data
taujk = {k:v for k,v in taujk.iteritems() if v < 50}
taujl = {k:v for k,v in taujl.iteritems() if v < 200}

# ------------------- Variables ----------------------
# Total annual cost (10^6 $)
TAC = mod.addVar(name="TAC")
# Total cost of feedstock during period t (10^6 $)
Cfeed = mod.addVars(times,name="Cfeed")
# Total cost of transportation during period t (10^6 $)
Ctrans = mod.addVars(times,name="Ctrans")
# Total cost of inventory during period t (10^6 $)
Cinv = mod.addVars(times,name="Cinv")
# Total production cost during period t (10^6 $)
Cprod = mod.addVars(times,name="Cprod")
# Total annual capital cost (10^6 $ per year)
Ccapex = mod.addVar(name="Ccapex")

# Sale of product from biorefinery during period t (10^6 liters)
R = mod.addVars(prods,brfs,times,name="R")
# Amount of compound harvested during period t (10^6 tons)
H = mod.addVars(feeds,hs,times,name="H")

# Amount of compound sent along arc during period t (10^6 tons)
Fjl = mod.addVars(feeds,hs,brfs,times,name="Fjl")
Fjk = mod.addVars(feeds,hs,dpts,times,name="Fjk")
Fkl = mod.addVars(itrmds,dpts,brfs,times,name="Fkl")

# Consumption of compound at deppot/biorefinery during period (10^6 tons)
Gl_f = mod.addVars(feeds,brfs,brfs_techs,times,name="Gl_f")
Gl_i = mod.addVars(itrmds,brfs,brfs_techs,times,name="Gl_i")
Gk = mod.addVars(feeds,brfs,dpts_techs,times,name="Gk")

# Production of compound at biorefinery/depot using technology m during period (10^6 tons)
Pl = mod.addVars(prods,brfs,brfs_techs,times,name="Pl")
Pk = mod.addVars(itrmds,dpts,dpts_techs,times,name="Pk")

# Daily operating capacity (tons per day)
Ql = mod.addVars(brfs,brfs_techs,name="Ql")
Qk = mod.addVars(dpts,dpts_techs,name="Qk")

# Inventory of biorefinery at the end of period t (10^6 tons)
Sl = mod.addVars(prods,brfs,times,name="Sl")
Sk = mod.addVars(itrmds,dpts,times,name="Sk")
Sj = mod.addVars(feeds,hs,times,name="Sj")

# Binary selection - will equal 1 if technology m is selected at refinery l
Uk = mod.addVars(dpts,dpts_techs,vtype=GRB.BINARY,name="Uk")

mod.update()

# -------------- Equations -----------------------------------------

mod.setObjective(TAC,GRB.MINIMIZE)

# Feedstock cost (10^6 $)
for t in times:
    mod.addConstrs(quicksum(H[i,j,t]*llambda[i] for i in feeds for j in hs) == Cfeed[t] for t in times)
# Inventory cost (10^6 $)
for t in times:
    mod.addConstr((quicksum(Sj[i,j,t]*iota[i] for i in feeds for j in hs)+quicksum(Sk[i,k,t]*iota[i] for i in itrmds for k in dpts)) == Cinv[t])
# Production cost (10^6 $)    
for t in times:
    (mod.addConstr((quicksum(Gl_f[i,l,m,t]*mu[m] for i in feeds for l in brfs for m in brfs_techs) + quicksum(Gl_i[i,l,m,t]*mu[m]
                  for i in itrmds for m in brfs_techs) + quicksum(Gk[i,k,m,t]*mu[m] for i in feeds for k in dpts for m in dpts_tech))) == Cprod[t])  
# Transportation cost (10^6 $)    
for t in times:
    mod.addConstr(quicksum((kappaTf[i]+kappaTv[i]*taujl[(j,l)])*Fjl[i,j,l,t] for i in feeds for (j,l) in taujl for l in brfs \
                           for j in hs)+quicksum((kappaTf[i]+kappaTv[i]*taujk[(j,k)])*Fjk[i,j,k,t] for i in feeds \
                            for (j,k) in taujk for k in dpts for j in hs) + quicksum((kappaRf[i]+kappaRv[i]*taukl[(k,l)])*Fkl[i,k,l,t] \
                              for i in itrmds for (k,l) in taukl for k in dpts for l in brfs) == Ctrans[t])    
# Annualized capital cost based on daily plant capacity (10^6 $)
mod.addConstr(quicksum((zeta[m]*Ql[l,m]/2000) for m in brfs_techs for l in brfs)+quicksum((zeta[m]*Qk[k,m]/200) for m in dpts_techs for k in dpts) == Ccapex)

# Mass Balance Constraints
# Refinery inventories (10^6 liters)
for i in prods:
    for l in brfs:
        for t in times:
            if t > 1:
                mod.addConstr(Sl[i,l,t-1]*(1-gamma[i])+quicksum(P[i,l,m,t] for m in brfs_techs)-R[i,l,t] == Sl[i,l,t])
            elif t==1:
                mod.addConstr(Sl[i,l,t+3]*(1-gamma[i])+quicksum(P[i,l,m,t] for m in brfs_techs)-R[i,l,t] == Sl[i,l,t])
                
# Depot site inventories (10^6 tons)
for i in itrmds:
    for k in dpts:
        for t in times:
            if t > 1:
                mod.addConstr(Sk[i,k,t-1]*(1-gamma[i])-quicksum(Fkl[i,k,l,t] for l in brfs)+quicksum(Pk[i,k,m,t] for m in dpts_techs) == Sk[i,k,t])
            elif t==1:
                mod.addConstr(Sk[i,k,t+3]*(1-gamma[i])-quicksum(Fkl[i,k,l,t] for l in brfs)+quicksum(Pk[i,k,m,t] for m in dpts_techs) == Sk[i,k,t])
