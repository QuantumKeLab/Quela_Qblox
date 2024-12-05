import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
import numpy as np
from qblox_instruments import Cluster
from Modularize.support.UserFriend import *
from scipy.optimize import minimize
from Modularize.support import QDmanager
from Modularize.m14_SingleShot import SS_executor, Qubit_state_single_shot
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, coupler_zctrl
import matplotlib.pyplot as plt

try:
    from Modularize.analysis.OneShotAna import a_OSdata_analPlot
    mode = "AS"
except:
    mode = "WeiEn"
    
def RO_optimization(QD_agent:QDmanager,cluster,Fctrl,target_q:str='q0',run:bool=True,shots:int=10000,plot:bool=True,exp_label:int=0,IF:float=250e6):
    print("Start to optimize readout")
    ## Initialize variables
    rof = []
    rol = []
    err = []

    ## Create QD_agent instance
    # QD_agent = QDmanager()
    qubit_info = QD_agent.quantum_device.get_element(target_q)
    old_rof = qubit_info.clock_freqs.readout()
    old_rol = qubit_info.measure.pulse_amp()
    
    # qubit_info.rxy.duration(6e-8)
    print(qubit_info.rxy.duration())
    ## Initial readout fidelity and error calculation
    # info= SS_executor(QD_agent,cluster, Fctrl,target_q)
    info=SS_executor(QD_agent,cluster,Fctrl,target_q=target_q,execution=run,shots=shots,roAmp_modifier=1,plot=plot if repeat ==1 else False,exp_label=exp_label,IF= IF)#roAmp_modifier=ro_amp_scaling,
          
    RO_fidelity = info[2]
    print(RO_fidelity)
    err_init = 1 - RO_fidelity
    rof.append(old_rof)
    rol.append(old_rol)
    err.append(err_init)

    ## Define the goal function to minimize the error
    def goal_function(params):
        rof, rol = params
        qubit_info.clock_freqs.readout(rof) 
        qubit_info.measure.pulse_amp(rol) 
        info= SS_executor(QD_agent)
        RO_fidelity = info[2]
        return 1 - RO_fidelity

    ## Define a callback function to track optimization progress
    optimization_trace = {'rof': [], 'rol': [], 'infidelity': []}

    def callback(params):
        f_val = goal_function(params)
        optimization_trace['infidelity'].append(f_val)
        optimization_trace['rof'].append(params[0])
        optimization_trace['rol'].append(params[1])

    ## Set initial values for rof and rol
    init_value = [old_rof, old_rol]
    
    ## Define bounds for rof and rol
    bounds = [(0.1, 0.9),  # rof's range [0.1, 0.9]
          (0.1, 0.9)]  # rol's range [0.1, 0.9]
    
    ## Run the optimization using Nelder-Mead
    result = minimize(goal_function, init_value, method='L-BFGS-B', callback=callback, bounds=bounds,tol=2e-3,
                    options={'maxiter': 200, 'disp': True, 'return_all': True})

    ## Extract the results
    optimized_rof, optimized_rol = result.x
    rof.extend(optimization_trace['rof'])
    rol.extend(optimization_trace['rol'])
    err.extend(optimization_trace['infidelity'])

    ## Plot the error progression
    plt.plot(range(len(err)), err)
    plt.xlabel('Iteration')
    plt.ylabel('Error')
    plt.title('Error Optimization Progress')
    plt.show()

def RO_opti_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,target_q:str,run:bool=True,shots:int=10000,plot:bool=True,exp_idx:int=0,IF:float=250e6):#,ro_amp_adj:float=1):
      
    if run:
        RO_optimization(QD_agent,cluster,Fctrl,target_q=target_q,run=run,shots=shots,plot=plot,exp_label=exp_idx, IF=IF)
    cluster.reset()
    
    
    
if __name__ == "__main__":
    
    """ Fill in """
    execution:bool= True
    repeat:int = 1
    DRandIP = {"dr":"drke","last_ip":"242"}
    # ro_elements = ['q0']  
    ro_elements = {'q0':{"roAmp_factor":1}}     
    couplers = []
    
    """ Optional paras (don't use is better) """
    ro_atte_degrade_dB:int = 0 # multiple of 2 
    shot_num:int = 10000
    xy_IF = 250e6


    for qubit in ro_elements:
        for i in range(repeat):
            """ Preparations """
            slightly_print(f"The {i}th OS:")
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"]) #r"C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_9_30\DRKE#242_SumInfo.pkl"
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
            QD_agent.Notewriter.modify_DigiAtte_For(-ro_atte_degrade_dB, qubit, 'ro')

            """ Running """
            
            Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
            
            RO_opti_executor(QD_agent,cluster,Fctrl,target_q=qubit,run=execution,shots=shot_num,plot=True,exp_idx=i,IF=xy_IF)
           
            
           
            # if ro_elements[qubit]["ro_amp_factor"] !=1:
            #     keep = mark_input(f"Keep this RO amp for {qubit}?[y/n]")
            # else:
            #     keep = 'y'


            """ Storing """
            # if execution and repeat == 1:
            #     if keep.lower() in ['y', 'yes']:
            #         QD_agent.QD_keeper() 
            
            
            """ Close """
            print("Readout optimization done!")
            shut_down(cluster,Fctrl,Cctrl)