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
    