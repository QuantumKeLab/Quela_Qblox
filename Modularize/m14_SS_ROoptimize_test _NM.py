###from DJ
from Modularize.support import QDmanager
from Modularize.m14_SingleShot import SS_executor
import matplotlib.pyplot as plt

rof=[]
rol=[]
err=[]

QD_agent=QDmanager()
old_rof=QD_agent.ROF
old_rol=QD_agent.ROL

RO_fidelity=SS_executor(QD_agent)
err_init=1-RO_fidelity
rof.append(old_rof)
rol.append(old_rol)
err.append(err_init)

def NM(rof,rol,err):
    pass

cond:bool=True
while cond:
    new_rof, new_rol=NM(old_rof,old_rol,2*(err_init))
    QD_agent.ROF =new_rof
    QD_agent.ROL =new_rol
    RO_fidelity=SS_executor(QD_agent)
    old_rof=new_rof
    old_rol=new_rol
    err_a=1-RO_fidelity
    rof.append(old_rof)
    rol.appemd(old_rol)
    err.append(err_a)
    
    if err_a<0.05:
        cond=False

plt.plot(x,err)
plt.show()


###from JY
from scipy.optimize import minimize

optimization_trace={'infidelity':[]}
for i in range(len(init_key)):
    optimization_trace[init_key[i]]=[]
    
def callback(init_value):
    f_val=goal_function(init_value)
    optimization_trace['infidelity'].append(f_val)
    for i in range(len(init_key)):
        optimization_trace[init_key[i]].append(init_value[i])
        
result=minimize(goal_function, init_value, method='Nelder-Mead',callback=callback, tol=2e-3,
                options={'maxdev':200,
                         'disp':True,
                         'return_all':True},
                )
    
###Combine by GPT
import numpy as np
from scipy.optimize import minimize
from Modularize.support import QDmanager
from Modularize.m14_SingleShot import SS_executor
import matplotlib.pyplot as plt

# Initialize variables
rof = []
rol = []
err = []

# Create QD_agent instance
QD_agent = QDmanager()
old_rof = QD_agent.ROF
old_rol = QD_agent.ROL

# Initial readout fidelity and error calculation
RO_fidelity = SS_executor(QD_agent)
err_init = 1 - RO_fidelity
rof.append(old_rof)
rol.append(old_rol)
err.append(err_init)

# Define the goal function to minimize the error
def goal_function(params):
    rof, rol = params
    QD_agent.ROF = rof
    QD_agent.ROL = rol
    RO_fidelity = SS_executor(QD_agent)
    return 1 - RO_fidelity

# Define a callback function to track optimization progress
optimization_trace = {'rof': [], 'rol': [], 'infidelity': []}

def callback(params):
    f_val = goal_function(params)
    optimization_trace['infidelity'].append(f_val)
    optimization_trace['rof'].append(params[0])
    optimization_trace['rol'].append(params[1])

# Set initial values for rof and rol
init_value = [old_rof, old_rol]

# Run the optimization using Nelder-Mead
result = minimize(goal_function, init_value, method='Nelder-Mead', callback=callback, tol=2e-3,
                  options={'maxiter': 200, 'disp': True, 'return_all': True})

# Extract the results
optimized_rof, optimized_rol = result.x
rof.extend(optimization_trace['rof'])
rol.extend(optimization_trace['rol'])
err.extend(optimization_trace['infidelity'])

# Plot the error progression
plt.plot(range(len(err)), err)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Error Optimization Progress')
plt.show()
