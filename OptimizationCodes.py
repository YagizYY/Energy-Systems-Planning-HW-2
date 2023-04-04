#%%
import os

# %%
import pandas as pd
import numpy as np
from pulp import *
import openpyxl
import matplotlib.pyplot as plt
# %%
import sys
print(sys.prefix)

#%%
inflow = pd.read_excel("C:/Users/yagiz/Desktop/4-2/IE-453/Energy-Systems-Planning-HW2/HybridSystemData.xlsx", sheet_name="StreamFlow")

solar_radiation = pd.read_excel("C:/Users/yagiz/Desktop/4-2/IE-453/Energy-Systems-Planning-HW2/HybridSystemData.xlsx", sheet_name="SolarRadiation")

demand = pd.read_excel("C:/Users/yagiz/Desktop/4-2/IE-453/Energy-Systems-Planning-HW2/HybridSystemData.xlsx", sheet_name="Demand")

# some manipulations on datasets regarding units
inflow["m3"] = inflow["km^3"]*(1000**3)
solar_radiation["KWh/m2"] = solar_radiation["GWh/km^2"]
demand["KWh"] = (10**6)*demand["GWh"]

inflow= inflow.drop("km^3", axis=1)
solar_radiation=solar_radiation.drop("GWh/km^2", axis=1)
demand=demand.drop("GWh", axis=1)

#%%
############################# PART A #############################
# initialize model
model = LpProblem(name="MinimizeCost", sense=LpMinimize)

# Define decision variables
purchase_grid = LpVariable.dicts("purchase_grid", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S = LpVariable.dicts("S", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

release = LpVariable.dicts("release", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

hydro_used = LpVariable.dicts("hydro_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_used = LpVariable.dicts("solar_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_curtailed = LpVariable.dicts("solar_curtailed", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill = LpVariable.dicts("spill", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

#parameters -> conver to meters & kwh
solar_farm_cap = 8*(10**7) # m^2
s_max = 10**8 # hydroenergy water capacity m^3
g_max = 1.2*(10**6) # generator max capacity -> kw
alpha = .12 # efficiency of solar
gamma = .9 # efficiency of generator
transmission_lost = .05 # transmission lines energy lost
g = 9.8 # m/s^2
h = 100 # m
d = 1000 # kg/m^3
#price = [0.3 for i in range(56)]
price = 0.3 # cent/kwh

# Define Objective Function
model += lpSum([price*purchase_grid[i] for i in demand.index])

# Define Constraints

# (1)
for i in range(1, (len(demand))):
    model += S[i] == S[i-1] + inflow["m3"][i] - release[i] - spill[i]

# (2)
model += S[0] == s_max/2 + inflow["m3"][0] - release[0] - spill[0]

# (3)
model += S[55] == s_max/2

# (4)
for i in demand.index:
    model += hydro_used[i] == release[i]*g*h*d*gamma/(3600*(10**3))

# (5)
for i in demand.index:
    model += solar_curtailed[i]+solar_used[i] == solar_radiation["KWh/m2"][i]*solar_farm_cap*alpha

# (6)
for i in demand.index:
    model += solar_used[i] + hydro_used[i]*(1-transmission_lost)+purchase_grid[i] >= demand["KWh"][i] 

# (7)
for i in demand.index:
    model += hydro_used[i] <= g_max*3 

# (8)
for i in demand.index:
    model += S[i] <= s_max


# Solve model
model.solve()


with open("decision_variable_values.txt", "w") as f:
    for v in model.variables():
        f.write(f"{v.name} = {v.varValue}\n")

print("Objective function value:", pulp.value(model.objective))

# Create a list of variable names and sort them numerically
var_names = sorted(
    [v.name for v in model.variables()],
    key=lambda x: (int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else -1, x),
)

# Create an empty dictionary to hold the values of the decision variables
var_dict = {}

# Loop through the decision variables and populate the dictionary
for v in model.variables():
    var_dict[v.name] = v.varValue

# Convert the dictionary to a pandas dataframe and sort the rows by variable name
df = pd.DataFrame.from_dict(var_dict, orient="index", columns=["Value"]).loc[var_names]

with open('decision_variable_outputs.txt', 'w') as f:
    for v in model.variables():
        f.write(f"{v.name}: {v.varValue}\n")


#%%
############################# PART C-1  #############################
##### NO RESERVOIR #####
# initialize model
model = LpProblem(name="MinimizeCost", sense=LpMinimize)

# Define decision variables
purchase_grid = LpVariable.dicts("purchase_grid", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

release = LpVariable.dicts("release", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

hydro_used = LpVariable.dicts("hydro_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_used = LpVariable.dicts("solar_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_curtailed = LpVariable.dicts("solar_curtailed", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill = LpVariable.dicts("spill", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

#parameters -> conver to meters & kwh
solar_farm_cap = 8*(10**7) # m^2
s_max = 10**8 # hydroenergy water capacity m^3
g_max = 1.2*(10**6) # generator max capacity -> kw
alpha = .12 # efficiency of solar
gamma = .9 # efficiency of generator
transmission_lost = .05 # transmission lines energy lost
g = 9.8 # m/s^2
h = 100 # m
d = 1000 # kg/m^3
#price = [0.3 for i in range(56)]
price = 0.3 # cent/kwh

# Define Objective Function
model += lpSum([price*purchase_grid[i] for i in demand.index])

# Define Constraints

# (1)
for i in demand.index:
    model += release[i] + spill[i] == inflow["m3"][i] 

# (2)
for i in demand.index:
    model += hydro_used[i] == release[i]*g*h*d*gamma/(3600*(10**3))

# (4)
for i in demand.index:
    model += solar_curtailed[i]+solar_used[i] == solar_radiation["KWh/m2"][i]*solar_farm_cap*alpha

# (5)
for i in demand.index:
    model += solar_used[i] + hydro_used[i]*(1-transmission_lost)+purchase_grid[i] >= demand["KWh"][i] 

# (6)
for i in demand.index:
    model += hydro_used[i] <= g_max*3 

# Solve model
model.solve()


with open("decision_variable_values.txt", "w") as f:
    for v in model.variables():
        f.write(f"{v.name} = {v.varValue}\n")

print("Objective function value:", pulp.value(model.objective))

# Create a list of variable names and sort them numerically
var_names = sorted(
    [v.name for v in model.variables()],
    key=lambda x: (int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else -1, x),
)

# Create an empty dictionary to hold the values of the decision variables
var_dict = {}

# Loop through the decision variables and populate the dictionary
for v in model.variables():
    var_dict[v.name] = v.varValue

# Convert the dictionary to a pandas dataframe and sort the rows by variable name
df = pd.DataFrame.from_dict(var_dict, orient="index", columns=["Value"]).loc[var_names]

with open('decision_variable_outputs.txt', 'w') as f:
    for v in model.variables():
        f.write(f"{v.name}: {v.varValue}\n")
 
# %%
############################# PART C-2 #############################
# initialize model
model = LpProblem(name="MinimizeCost", sense=LpMinimize)

# Define decision variables
purchase_grid = LpVariable.dicts("purchase_grid", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S = LpVariable.dicts("S", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

release = LpVariable.dicts("release", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

hydro_used = LpVariable.dicts("hydro_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_used = LpVariable.dicts("solar_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_curtailed = LpVariable.dicts("solar_curtailed", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill = LpVariable.dicts("spill", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

#parameters -> conver to meters & kwh
solar_farm_cap = 8*(10**7) # m^2
s_max = 10**8 # hydroenergy water capacity m^3
g_max = 0.015*(10**6) # generator max capacity -> kw
alpha = .12 # efficiency of solar
gamma = .9 # efficiency of generator
transmission_lost = .05 # transmission lines energy lost
g = 9.8 # m/s^2
h = 100 # m
d = 1000 # kg/m^3
#price = [0.3 for i in range(56)]
price = 0.3 # cent/kwh

# Define Objective Function
model += lpSum([price*purchase_grid[i] for i in demand.index])

# Define Constraints

# (1)
for i in range(1, (len(demand))):
    model += S[i] == S[i-1] + inflow["m3"][i] - release[i] - spill[i]

# (2)
model += S[0] == s_max/2 + inflow["m3"][0] - release[0] - spill[0]

# (3)
model += S[55] == s_max/2

# (4)
for i in demand.index:
    model += hydro_used[i] == release[i]*g*h*d*gamma/(3600*(10**3))

# (5)
for i in demand.index:
    model += solar_curtailed[i]+solar_used[i] == solar_radiation["KWh/m2"][i]*solar_farm_cap*alpha

# (6)
for i in demand.index:
    model += solar_used[i] + hydro_used[i]*(1-transmission_lost)+purchase_grid[i] >= demand["KWh"][i] 

# (7)
for i in demand.index:
    model += hydro_used[i] <= g_max*3 

# (8)
for i in demand.index:
    model += S[i] <= s_max


# Solve model
model.solve()


with open("decision_variable_values.txt", "w") as f:
    for v in model.variables():
        f.write(f"{v.name} = {v.varValue}\n")

print("Objective function value:", pulp.value(model.objective))

# Create a list of variable names and sort them numerically
var_names = sorted(
    [v.name for v in model.variables()],
    key=lambda x: (int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else -1, x),
)

# Create an empty dictionary to hold the values of the decision variables
var_dict = {}

# Loop through the decision variables and populate the dictionary
for v in model.variables():
    var_dict[v.name] = v.varValue

# Convert the dictionary to a pandas dataframe and sort the rows by variable name
df = pd.DataFrame.from_dict(var_dict, orient="index", columns=["Value"]).loc[var_names]

with open('decision_variable_outputs.txt', 'w') as f:
    for v in model.variables():
        f.write(f"{v.name}: {v.varValue}\n")


# %%
############################# PART D #############################
########### Calculate Capacity Factor by using the model in A #######

actual_output_solar = df[df.index.str.startswith('solar_curtailed')].sum() + df[df.index.str.startswith('solar_used')].sum()

potential_output_solar = solar_radiation["KWh/m2"].max()*alpha*solar_farm_cap*56

capacity_factor_solar = actual_output_solar / potential_output_solar


actual_output_hydro = df[df.index.str.startswith('hydro_used')].sum()

potential_output_hydro = min(s_max*d*g*h*gamma/(3600*(10**3)), g_max*3*gamma)*56

capacity_factor_hydro = actual_output_hydro / potential_output_hydro

#%%
############################# PART E #############################
# initialize model
model = LpProblem(name="MinimizeCost", sense=LpMinimize)

# Define decision variables
purchase_grid = LpVariable.dicts("purchase_grid", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S_U = LpVariable.dicts("S_U", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S_L = LpVariable.dicts("S_L", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

release = LpVariable.dicts("release", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

pump = LpVariable.dicts("pump", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

hydro_used = LpVariable.dicts("hydro_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_used = LpVariable.dicts("solar_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_pumped = LpVariable.dicts("solar_pumped", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_curtailed = LpVariable.dicts("solar_curtailed", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill_L = LpVariable.dicts("spill_L", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill_U = LpVariable.dicts("spill_U", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

#parameters -> conver to meters & kwh
solar_farm_cap = 8*(10**7) # m^2
s_max_u = .1*(10**9) # hydroenergy water of upper reservoir capacity m^3
s_max_l = .05*(10**9) # hydroenergy water of lower reservoir capacity m3
g_max = 1.2*(10**6) # generator max capacity -> kw
alpha = .12 # efficiency of solar
gamma = .9 # efficiency of generator
transmission_lost = .05 # transmission lines energy lost
g = 9.8 # m/s^2
h = 100 # m
d = 1000 # kg/m^3
#price = [0.3 for i in range(56)]
price = 0.3 # cent/kwh

# Define Objective Function
model += lpSum([price*purchase_grid[i] for i in demand.index])

# Define Constraints

# (1)
for i in range(1, (len(demand))):
    model += S_U[i] == S_U[i-1] + inflow["m3"][i] - release[i] - spill_U[i] + pump[i]

# (2)
model += S_U[0] == s_max_u/2 + inflow["m3"][0] - release[0] - spill_U[0] + pump[0]

# (3)
for i in range(1, (len(demand))):
    model += S_L[i] == S_L[i-1] + release[i] - spill_L[i] - pump[i]

# (4)
model += S_L[0] == 0 + release[0] - spill_L[0] - pump[0]

# (5)
for i in demand.index:
    model += solar_curtailed[i]+solar_used[i] + solar_pumped[i] == solar_radiation["KWh/m2"][i]*solar_farm_cap*alpha

# (6)
for i in demand.index:
    model += hydro_used[i] == release[i]*g*h*d*gamma/(3600*(10**3))

# (7)
for i in demand.index:
    model += hydro_used[i] <= g_max*3 

# (8)
for i in demand.index:
    model += pump[i]*d*g*h/(gamma*(3600*(10**3))) == solar_pumped[i]

# (9)
for i in demand.index:
    model += solar_used[i] + hydro_used[i]*(1-transmission_lost)+purchase_grid[i] >= demand["KWh"][i] 


# (10)
for i in demand.index:
    model += S_U[i] <= s_max_u

# (11)
for i in demand.index:
    model += S_L[i] <= s_max_l


# Solve model
model.solve()

with open("decision_variable_values.txt", "w") as f:
    for v in model.variables():
        f.write(f"{v.name} = {v.varValue}\n")

print("Objective function value:", pulp.value(model.objective))

# Create a list of variable names and sort them numerically
var_names = sorted(
    [v.name for v in model.variables()],
    key=lambda x: (int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else -1, x),
)

# Create an empty dictionary to hold the values of the decision variables
var_dict = {}

# Loop through the decision variables and populate the dictionary
for v in model.variables():
    var_dict[v.name] = v.varValue

# Convert the dictionary to a pandas dataframe and sort the rows by variable name
df = pd.DataFrame.from_dict(var_dict, orient="index", columns=["Value"]).loc[var_names]

with open('decision_variable_outputs.txt', 'w') as f:
    for v in model.variables():
        f.write(f"{v.name}: {v.varValue}\n")

#%%
############################# PART F #############################
# initialize model
model = LpProblem(name="MinimizeCost", sense=LpMinimize)

# Define decision variables
purchase_grid = LpVariable.dicts("purchase_grid", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S_U = LpVariable.dicts("S_U", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

S_L = LpVariable.dicts("S_L", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

release = LpVariable.dicts("release", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

pump = LpVariable.dicts("pump", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

hydro_used = LpVariable.dicts("hydro_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_used = LpVariable.dicts("solar_used", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_pumped = LpVariable.dicts("solar_pumped", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

solar_curtailed = LpVariable.dicts("solar_curtailed", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill_L = LpVariable.dicts("spill_L", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

spill_U = LpVariable.dicts("spill_U", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

shift = LpVariable.dicts("shift", [i for i in demand.index],lowBound=0, upBound=None, cat="Continuous")

#parameters -> conver to meters & kwh
solar_farm_cap = 8*(10**7) # m^2
s_max_u = .1*(10**9) # hydroenergy water of upper reservoir capacity m^3
s_max_l = .05*(10**9) # hydroenergy water of lower reservoir capacity m3
g_max = 1.2*(10**6) # generator max capacity -> kw
alpha = .12 # efficiency of solar
gamma = .9 # efficiency of generator
transmission_lost = .05 # transmission lines energy lost
g = 9.8 # m/s^2
h = 100 # m
d = 1000 # kg/m^3
#price = [0.3 for i in range(56)]
price = 0.3 # cent/kwh

# Define Objective Function
model += lpSum([price*purchase_grid[i] for i in demand.index])

# Define Constraints

# (1)
for i in range(1, (len(demand))):
    model += S_U[i] == S_U[i-1] + inflow["m3"][i] - release[i] - spill_U[i] + pump[i]

# (2)
model += S_U[0] == s_max_u/2 + inflow["m3"][0] - release[0] - spill_U[0] + pump[0]

# (3)
for i in range(1, (len(demand))):
    model += S_L[i] == S_L[i-1] + release[i] - spill_L[i] - pump[i]

# (4)
model += S_L[0] == 0 + release[0] - spill_L[0] - pump[0]

# (5)
for i in demand.index:
    model += solar_curtailed[i]+solar_used[i] + solar_pumped[i] == solar_radiation["KWh/m2"][i]*solar_farm_cap*alpha

# (6)
for i in demand.index:
    model += hydro_used[i] == release[i]*g*h*d*gamma/(3600*(10**3))

# (7)
for i in demand.index:
    model += hydro_used[i] <= g_max*3 

# (8)
for i in demand.index:
    model += pump[i]*d*g*h/(gamma*(3600*(10**3))) == solar_pumped[i]

# (9)
model += solar_used[0] + hydro_used[0]*(1-transmission_lost)+purchase_grid[0] >= demand["KWh"][0] - shift[0] 

# (10)
for i in range(1, (len(demand))):
    model += solar_used[i] + hydro_used[i]*(1-transmission_lost)+purchase_grid[i] >= demand["KWh"][i] - shift[i] + shift[i-1] 

# (11)
for i in demand.index:
    model += shift[i] <= demand["KWh"][i]

# (12)
model += shift[55] == 0

# (13)
for i in demand.index:
    model += S_U[i] <= s_max_u

# (14)
for i in demand.index:
    model += S_L[i] <= s_max_l


# Solve model
model.solve()

with open("decision_variable_values.txt", "w") as f:
    for v in model.variables():
        f.write(f"{v.name} = {v.varValue}\n")

print("Objective function value:", pulp.value(model.objective))

# Create a list of variable names and sort them numerically
var_names = sorted(
    [v.name for v in model.variables()],
    key=lambda x: (int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else -1, x),
)

# Create an empty dictionary to hold the values of the decision variables
var_dict = {}

# Loop through the decision variables and populate the dictionary
for v in model.variables():
    var_dict[v.name] = v.varValue

# Convert the dictionary to a pandas dataframe and sort the rows by variable name
df = pd.DataFrame.from_dict(var_dict, orient="index", columns=["Value"]).loc[var_names]

with open('decision_variable_outputs.txt', 'w') as f:
    for v in model.variables():
        f.write(f"{v.name}: {v.varValue}\n")
