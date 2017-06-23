# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 14:46:39 2017

@author: boigenzahn
"""
from gurobipy import *
import pandas as pd
import numpy as np

file = pd.ExcelFile('..\gamsdir\projdir\Task1b.xlsx')

m = Model('Task1b_a')

# Parse once and split the dataframe, not the Excel Sheet
# set_index as argument of parse - also from unicode argument
# Harvest sites
hs_df = file.parse(header=1,skipfooter=677,parse_cols='M')
u_hs = hs_df.values.flatten()
hs = [str(x) for x in u_hs]
#print(hs)    

#PW availability in 10^3 tons - Using PW3
alpha_df = file.parse(header=1, skipfooter=677,parse_cols='B,E', converters={'*':float})
u_alpha = alpha_df.set_index('Sites').T.to_dict('list')
alpha = {str(k): v for k, v in u_alpha.items()}
#print(alpha)

#Distance to biorefinery 1 in miles
tau_df = file.parse(header=2,skip_footer=677,parse_cols='M:N')
u_tau = tau_df.set_index('Sites').T.to_dict('list')
tau = {str(k): v for k, v in u_tau.items()}
#print(tau)

#Constants
d = 47.25  #Ethanol demand in 10^6 liters
y = 283.906     #wood to EtOH yield in liters/ton
ft = 4.29   #Fixed transportation cost in $ per ton
vt = 0.0737 #Variable transportaiton cost in $ per ton mile

# Unit cost of travel from site j ($/ton)
gamma = {}
for k, v in tau.iteritems():
    gamma[k] = vt*v[0] + ft
#print(gamma)

#amount shipped from harvest site j
F = m.addVars(hs,name="quantity")

# Transportation cost variable
TC = m.addVar(name="TC")

m.update()

m.setObjective(TC, GRB.MINIMIZE)

#cost..        TC =e= sum(j, gamma(j)*F(j)) ;    
m.addConstr(quicksum(gamma[j]*F[j] for j in hs) == TC)

#supply(j)..     F(j) =l= alpha(j) ;
m.addConstrs(F[j] <= alpha[j][0] for j in hs)

#demand..        sum(j, F(j)) =g= (d*1000)/y ;  
m.addConstr(quicksum(F[j] for j in hs) >= (d*1000)/y)

m.optimize()