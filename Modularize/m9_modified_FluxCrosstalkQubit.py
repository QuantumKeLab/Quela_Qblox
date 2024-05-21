import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
import matplotlib.pyplot as plt
from numpy import array, linspace
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from utils.tutorial_analysis_classes import QubitFluxSpectroscopyAnalysis
from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset
from Modularize.support.ReadResults import plot_QbFlux, Zline_Crosstalk
from Modularize.support.Pulse_schedule_library import Zline_crosstalk_sche, set_LO_frequency, pulse_preview
import numpy as np
import xarray as xr

def update_2Ddict(dict, key_a, key_b, val):
    if key_a in dict:
        dict[key_a].update({key_b:val})
    else:
        dict.update({key_a:{key_b:val}})

def ZlineCrosstalk_two_tone_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,Z_amp_start:float,Z_amp_end:float,IF:int=100e6,n_avg:int=1000,crosstalk_Z_points:int=20,target_Z_points:int=20,run:bool=True,q:str='q1',z:str='q1',Experi_info={},get_data_path:bool=False,factor:float=0.02,compensated_dict:dict={}):
    print("ZlineCrosstalk 2tone start")
    sche_func = Zline_crosstalk_sche
    qubit_info = QD_agent.quantum_device.get_element(q)
    original_f01 = qubit_info.clock_freqs.f01()
    print(original_f01)

    LO= qubit_info.clock_freqs.f01()+IF
    print(qubit_info.clock_freqs.f01()*1e-9)
    set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=LO)
    
    q_Z_bias = ManualParameter(name="q_Z", unit="V", label="q_Z bias")
    q_Z_bias.batched = False
    z_Z_bias = ManualParameter(name="z_Z", unit="V", label="z_Z bias")
    z_Z_bias.batched = False
    
    # temperature quard
    if Z_amp_end > 0.4:
        Z_amp_end = 0.4
    elif Z_amp_end < -0.4:
        Z_amp_end = -0.4
    else:
        pass

    if Z_amp_start > 0.4:
        Z_amp_start = 0.4
    elif Z_amp_start < -0.4:
        Z_amp_start = -0.4
    else:
        pass 

    q_Z_samples = linspace(factor*Z_amp_start,factor*Z_amp_end,target_Z_points)
    z_Z_samples = linspace(Z_amp_start,Z_amp_end,crosstalk_Z_points)

    sched_kwargs = dict(   
        q=q,
        z=z,

        # pi_amp={str(q):QD_agent.Notewriter.get_2tone_piampFor(q)},
        # pi_dura=48e-6,

        pi_amp={str(q):0.5*qubit_info.rxy.amp180()},
        pi_dura=2*qubit_info.rxy.duration(),

        q_Z_amp=q_Z_bias,
        z_Z_amp=z_Z_bias,
        R_amp={str(q):qubit_info.measure.pulse_amp()},
        R_duration={str(q):qubit_info.measure.pulse_duration()},
        R_integration={str(q):qubit_info.measure.integration_time()},
        R_inte_delay=qubit_info.measure.acq_delay(),
        compensation=compensated_dict
    )
    exp_kwargs= dict(q_Z_amp=['start '+'%E' %q_Z_samples[0],'end '+'%E' %q_Z_samples[-1]],
                     z_Z_amp=['start '+'%E' %z_Z_samples[0],'end '+'%E' %z_Z_samples[-1]],)
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=sched_kwargs,
            real_imag=False,
            batched=False,
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([q_Z_bias,z_Z_bias])
        meas_ctrl.setpoints_grid((q_Z_samples,z_Z_samples))
        qs_ds = meas_ctrl.run("ZlineCrosstalk-two-tone")
        print(qs_ds)
        
        # Save the raw data into netCDF
        if get_data_path:
            path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='Z2tone',get_data_loc=get_data_path)
        else:
            path = ''
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=qs_ds,qb=q,exp_type='Z2tone',get_data_loc=get_data_path)
                
        show_args(exp_kwargs, title="ZlineCrosstalk_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        
    else:
        print('not running')
        sweep_para1= array(q_Z_samples[:2])
        sweep_para2= array(z_Z_samples[:2])
        sched_kwargs['q_Z_amp']= sweep_para1.reshape(sweep_para1.shape or (1,))[1]
        sched_kwargs['z_Z_amp']= sweep_para2.reshape(sweep_para2.shape or (1,))[1]
        pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
        
        
        show_args(exp_kwargs, title="ZlineCrosstalk_two_tone_kwargs: Meas.qubit="+q)
        if Experi_info != {}:
            show_args(Experi_info(q))
        path = ''

    qubit_info.clock_freqs.f01(original_f01)

    return path


def fluxcrosstalkQubit_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,specific_qubits:str,specific_zlines:str,run:bool=True,z_shifter:float=0,crosstalk_zpts:int=20,target_zpts:int=20,z_span=0.4,crosstalk_factor:float=0.02,Compensation_dict:dict={}):
    center = QD_agent.Fluxmanager.get_sweetBiasFor(target_q=specific_qubits)
    if run:
        Fctrl[specific_qubits](center)
        Fctrl[specific_zlines](0.)
        nc_path = ZlineCrosstalk_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=-z_span+z_shifter,crosstalk_Z_points=crosstalk_zpts,target_Z__points=target_zpts,Z_amp_end=z_span+z_shifter,q=specific_qubits,z=specific_zlines,run=True,get_data_path=True,factor=crosstalk_factor,compensated_dict=Compensation_dict)
        crosstalk = Zline_Crosstalk(QD_agent,nc_path,specific_qubits)
        reset_offset(Fctrl)
        if not Compensation_dict:
            permission = input("Does it look good ? [y/n]") 
            if permission.lower() in ['y','yes']:
                print("##################",nc_path)
                return crosstalk
            else:
                print("bad")
                qubit = QD_agent.quantum_device.get_element(specific_qubits)
                check_again.append(qubit)
                return 0
        else:
            print("has been compensated")
    else:
        _= ZlineCrosstalk_two_tone_spec(QD_agent,meas_ctrl,Z_amp_start=-z_span+z_shifter,crosstalk_Z_points=crosstalk_zpts,target_Z__points=target_zpts,Z_amp_end=z_span+z_shifter,q=specific_qubits,z=specific_zlines,run=False,get_data_path=True,factor=crosstalk_factor,compensated_dict=Compensation_dict)
        return 0


if __name__ == "__main__":
    
    """ Fill in """
    execution = True
    DRandIP = {"dr":"dr3","last_ip":"13"}
    ro_elements = ['q0','q1',]
    Z_matrix = np.identity(len(ro_elements))
    z_shifter = 0 # V
    compensation_dict = {}

    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    if ro_elements == 'all':
        ro_elements = list(Fctrl.keys())

    
    """ Running """
    check_again =[]
    for i, qubit in enumerate(ro_elements):
        for j, zline in enumerate(ro_elements):
            if (qubit==zline):
                continue
            else:
                init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
                Z_matrix[i,j] = fluxcrosstalkQubit_executor(QD_agent,meas_ctrl,qubit,zline,run=execution,z_shifter=z_shifter,crosstalk_zpts = 20, target_zpts=20, z_span=0.4, crosstalk_factor=0.02, Compensation_dict=compensation_dict)
                cluster.reset()
    if not check_again:
        print("Z_matrix is")
        print(Z_matrix)
        print("\n")
        compensation_matrix = np.linalg.inv(Z_matrix)
        print("compensation_matrix is")
        print(compensation_matrix)
        print("\n")
        for i, qubit in enumerate(ro_elements):
            for j, zline in enumerate(ro_elements):
                update_2Ddict(compensation_dict, qubit, zline, compensation_matrix[j,i])
        for i, qubit in enumerate(ro_elements):
            for j, zline in enumerate(ro_elements):
                if (qubit==zline):
                    continue
                else:
                    init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
                    fluxcrosstalkQubit_executor(QD_agent,meas_ctrl,qubit,zline,run=execution,z_shifter=z_shifter,crosstalk_zpts = 20, target_zpts=20, z_span=0.4, crosstalk_factor=0.02, Compensation_dict=compensation_dict)
                    cluster.reset()
        permission = input("Do you want to save the result ? [y/n]") 
        if permission.lower() in ['y','yes']:
            QD_agent.Fluxmanager.巴拉巴拉(target_q=qubit,compensation=compensation_dict)
        else:
            print("Why not ?")


        
    else:
        print(f"qubits to check again: {check_again}")
        compensation_dict={}



    """ Storing """
    if  execution:
        QD_agent.QD_keeper()
    


    """ Close """
    print('Flux qubit done!')
    shut_down(cluster,Fctrl)
    


