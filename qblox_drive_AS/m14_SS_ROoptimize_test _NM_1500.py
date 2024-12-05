import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from numpy import array, linspace, median, std
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support.Pulse_schedule_library import Qubit_state_single_shot_plot
from qblox_drive_AS.support import QDmanager, Data_manager,init_system_atte, init_meas, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import Qubit_SS_sche, set_LO_frequency, pulse_preview, Qubit_state_single_shot_fit_analysis
from qblox_drive_AS.SingleReadout.m14_SingleShot import SS_executor
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

try:
    mode = "AS"
except:
    mode = "WeiEn"
    
def RO_optimization(QD_agent:QDmanager,cluster,Fctrl,target_q:str='q0',run:bool=True,shots:int=10000,plot:bool=True,exp_label:int=0,IF:float=250e6):
    print("Start to optimize readout")
   
    rof = []
    rol = []
    err = []

    qubit_info = QD_agent.quantum_device.get_element(target_q)
    old_rof = qubit_info.clock_freqs.readout()
    old_rol = qubit_info.measure.pulse_amp()
    
    qubit_info.rxy.duration(6e-8)
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


def RO_opti_executor(QD_agent:QDmanager,cluster:Cluster,Fctrl:dict,target_q:str,run:bool=True,shots:int=10000,plot:bool=True,exp_idx:int=0,IF:float=250e6):#,ro_amp_adj:float=1):
      
    if run:
        RO_optimization(QD_agent,cluster,Fctrl,target_q=target_q,run=run,shots=shots,plot=plot,exp_label=exp_idx, IF=IF)
    cluster.reset()
    
    
if __name__ == '__main__':
    

    """ Fill in """
    execute:bool = 1
    repeat:int = 1
    DRandIP = {"dr":"drke","last_ip":"242"}
    ro_elements = {'q0':{}}#"roAmp_factor":1
    couplers = []


    """ Optional paras (don't use is better) """
    ro_atte_degrade_dB:int = 0 # multiple of 2 
    shot_num:int = 10000
    xy_IF = 250e6



    """ Iteration """

    for qubit in ro_elements:
        for i in range(repeat):
            start_time = time.time()

            """ Preparation """
            slightly_print(f"The {i}th OS:")
            QD_path =find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path)#,mode='l'
            QD_agent.Notewriter.modify_DigiAtte_For(-ro_atte_degrade_dB, qubit, 'ro')


            """ Running """
            # Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
            init_system_atte(QD_agent.quantum_device,list([qubit]),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'))
    
            # info = SS_executor(QD_agent,cluster,Fctrl,qubit,execution=execute,shots=shot_num,roAmp_modifier=1,plot=True if repeat ==1 else False,exp_label=i,IF=xy_IF)#,data_folder=r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_25\SS_overnight"
            RO_opti_executor(QD_agent,cluster,Fctrl,target_q=qubit,run=execute,shots=shot_num,plot=True,exp_idx=i,IF=xy_IF)
                
            """ Close """    
            # shut_down(cluster,Fctrl,Cctrl)
            end_time = time.time()
            slightly_print(f"time cose: {round(end_time-start_time,1)} secs")

   

        
    #    ## Define the goal function to minimize the error
    # def goal_function(params):
    #     rof, rol = params
    #     qubit_info.clock_freqs.readout= rof
    #     qubit_info.measure.pulse_amp= rol
    #     info=SS_executor(QD_agent,cluster,Fctrl,target_q=target_q,execution=run,shots=shots,roAmp_modifier=1,plot=plot if repeat ==1 else False,exp_label=exp_label,IF= IF)#roAmp_modifier=ro_amp_scaling,
    #     RO_fidelity = info[2]
    #     return 1 - RO_fidelity

    # ## Define a callback function to track optimization progress
    # optimization_trace = {'rof': [], 'rol': [], 'infidelity': []}

    # def callback(params):
    #     f_val = goal_function(params)
    #     optimization_trace['infidelity'].append(f_val)
    #     optimization_trace['rof'].append(params[0])
    #     optimization_trace['rol'].append(params[1])

    # ## Set initial values for rof and rol
    # init_value = [old_rof, old_rol]
    
    # ## Define bounds for rof and rol
    # bounds = [(0, 0.9),  # rof's range [0.1, 0.9]
    #       (0.5, 2)]  # rol's range [0.1, 0.9]
    
    # ## Run the optimization using Nelder-Mead
    # result = minimize(goal_function, init_value, method='L-BFGS-B', bounds=bounds,tol=2e-3,
    #                 options={'maxiter': 200, 'disp': True})#callback=callback, 

    # ## Extract the results
    # optimized_rof, optimized_rol = result.x
    # rof.extend(optimization_trace['rof'])
    # rol.extend(optimization_trace['rol'])
    # err.extend(optimization_trace['infidelity'])

    # ## Plot the error progression
    # plt.plot(range(len(err)), err)
    # plt.xlabel('Iteration')
    # plt.ylabel('Error')
    # plt.title('Error Optimization Progress')
    # plt.show()
