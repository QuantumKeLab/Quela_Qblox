import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from numpy import NaN
from numpy import array, linspace, sqrt
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from Modularize.support.UserFriend import *
from Modularize.support import QDmanager, Data_manager, cds,compose_para_for_multiplexing
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas, init_system_atte, shut_down, reset_offset, coupler_zctrl
from Modularize.support.StarkShift  import plot_ROopti
from Modularize.support.Pulse_schedule_library import StarkShift_sche, set_LO_frequency, pulse_preview
import xarray as xr
import numpy as np
from xarray import Dataset



def Stark_shift_spec(QD_agent:QDmanager,meas_ctrl:MeasurementControl,ro_elements:str='q0',ro_p_max:float=0.5,IF:int=200e6,xyf:float=0e9,xyf_span_Hz:float=400e6,n_avg:int=1000,p_points:int=40,f_points:int=60,run:bool=True,get_data_path:bool=False,analysis:bool=True):#,q:str='q1'):#,Experi_info={},get_data_path:bool=False):
    print("Stark shift start")
    trustable = True
    sche_func=StarkShift_sche

    analysis_result = {}
    original_xyfs = {}
    
    
    
    for q in ro_elements:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.integration_time(1.2e-6)
        qubit_info.measure.pulse_duration(1.2e-6)
        print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
        print("Pulse duration ",qubit_info.measure.pulse_duration()*1e6, "µs")
        print("Reset time ", qubit_info.reset.duration()*1e6, "µs")
        
        original_xyfs = qubit_info.clock_freqs.f01()
        print("original_xyfs",original_xyfs)

        if xyf == 0:
            xyf_highest = original_xyfs+IF
        else:
            xyf_highest = xyf + IF
        qubit_info.clock_freqs.f01(NaN)
        set_LO_frequency(QD_agent.quantum_device,q=q,module_type='drive',LO_frequency=xyf_highest)
        print("LO_frequency",xyf_highest)
    
    f01_samples = linspace(xyf_highest-xyf_span_Hz,xyf_highest,f_points)#
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    ro_p_samples = linspace(0,ro_p_max,p_points)


    for q in ro_elements:
        spec_pulse_amp = QD_agent.Notewriter.get_2tone_piampFor(q)
    print('pi duration:',qubit_info.rxy.duration())
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp_2=ro_pulse_amp,
        pi_amp={str(q):spec_pulse_amp},
        pi_dura=qubit_info.rxy.duration(),
        R_amp=compose_para_for_multiplexing(QD_agent,ro_elements,'1'),
        R_duration=compose_para_for_multiplexing(QD_agent,ro_elements,'3'),
        R_duration_2={str(q):8e-6},
        R_integration=compose_para_for_multiplexing(QD_agent,ro_elements,'4'),
        R_inte_delay=compose_para_for_multiplexing(QD_agent,ro_elements,'2'),
        # correlate_delay:float=1200e-9, # 

    )
    
    if run:
        gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=True,
            batched=True,
            # num_channels=len(list(XYFs.keys())),
        )
        QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((f01_samples,ro_p_samples))
        qs_ds = meas_ctrl.run("Stark-Shift")
        
        dict_ = {}
        for idx, q in enumerate(ro_elements):#why use ro_elements if there is no following code?
            # print(q)
            # print(type(q),type(ro_elements))
            freq_values = 2*ro_p_samples.shape[0]* list(f01_samples)
            i_data = array(qs_ds[f'y{2*idx}']).reshape(ro_p_samples.shape[0],f01_samples.shape[0])
            q_data = array(qs_ds[f'y{2*idx+1}']).reshape(ro_p_samples.shape[0],f01_samples.shape[0])    
            
            dict_[q] = (["mixer","amp","freq"],array([i_data,q_data]))
            #dict_[f"{q}_freq"] = (["mixer", "amp", "freq"], array(freq_values).reshape(2, ro_p_samples.shape[0], f01_samples.shape[0]))
            
        rfs_ds = Dataset(dict_,coords={"mixer":array(["I","Q"]),"amp":ro_p_samples,"freq":f01_samples})
        rfs_ds.attrs["execution_time"] = Data_manager().get_time_now()
        
        # Save the raw data into netCDF
        if get_data_path:
            path = Data_manager().save_raw_data(QD_agent=QD_agent,ds=rfs_ds,qb=q,exp_type='starkshift',get_data_loc=get_data_path)
        else:
            path = ''
            Data_manager().save_raw_data(QD_agent=QD_agent,ds=rfs_ds,qb=q,exp_type='starkshift',get_data_loc=get_data_path)
        
        if analysis:
            try:
                analysis_result[q] = QubitFluxSpectroscopyAnalysis(tuid=qs_ds.attrs["tuid"], dataset=qs_ds).run()
            except:
                analysis_result[q] = {}
                print("Qb vs Flux fitting failed! Raw data had been saved.")
                trustable = False
   
    else: ##Multi readout
        n_s = 2  # 這表示要選擇每個元素的前 `n_s` 個樣本
        

        # 遍歷 `ro_elements`，這裡 `ro_elements` 是一個列表，包含量子比特名稱
        sweep_para2 = array([ro_p_samples[0],ro_p_samples[-1]])# 選取 ro_p_samples 的前 2 個數據作為掃描參數

        # 更新 `spec_sched_kwargs`，將頻率和掃描參數傳入
        spec_sched_kwargs['frequencies'] =  f01_samples[:n_s]
        spec_sched_kwargs['R_amp_2'] = ro_p_samples[int(ro_p_samples.shape[0]/2)]

        pulse_preview(QD_agent.quantum_device, sche_func, spec_sched_kwargs) # 預覽脈衝
       
        rfs_ds = "" # 如果不需要返回 rfs_ds，將其設為空字符串

    return path, trustable#rfs_ds


def update_by_StarkShift(QD_agent:QDmanager,correct_results:dict,target_q:str):
    """
    correct_results dict in the form: {"xyf":float,"sweet_bias":float}  #?
    """
    qubit = QD_agent.quantum_device.get_element(target_q)
    # qubit.clock_freqs.f01(correct_results["xyf"])
    qubit.measure.pulse_amp(correct_results[""])
    # QD_agent.Fluxmanager.check_offset_and_correctFor(target_q=target_q,new_offset=correct_results["sweet_bias"])
    # QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=target_q,bias=correct_results["sweet_bias"])

def StarkShift_executor(QD_agent:QDmanager,meas_ctrl:MeasurementControl,specific_qubits:str,run:bool=True,p_pts:int=20,max_p:float=0.15,fpts:int=40,f_sapn_Hz=400e6,avg_times:int=1000,xy_IF:float=200e6):
    if run:
        
        # Fctrl[specific_qubits](center)
        nc_path, trustable= Stark_shift_spec(QD_agent,meas_ctrl,ro_p_max=max_p,p_points=p_pts,f_points=fpts,ro_elements=specific_qubits,run=run,xyf_span_Hz=f_sapn_Hz,IF=xy_IF,n_avg=avg_times,get_data_path=True)
        # delta=1
        # g=1
        # n_critic_origin=delta**2/(4*g**2)
        Kai_eff=1
        AC_Shift=1
        n_critic=AC_Shift/(2* Kai_eff) 
        
        reset_offset(Fctrl)
        if trustable:
            plot_ROopti(QD_agent,nc_path)
            permission = mark_input("Update the QD with this result ? [y/n]") 
            if permission.lower() in ['y','yes']:
                return trustable, {"xyf":results[specific_qubits].quantities_of_interest["freq_0"].nominal_value,"sweet_bias":results[specific_qubits].quantities_of_interest["offset_0"].nominal_value+center}
            else:
                return False, {}
        else:
            plot_ROopti(QD_agent,nc_path)
            trustable = False
            return False, {}

    else:
        path, trustable= Stark_shift_spec(QD_agent,meas_ctrl,ro_elements=ro_elements,ro_p_max=ro_p_max,p_points=p_pts,f_points=fpts,run=run,get_data_path=True,xyf_span_Hz=f_sapn_Hz,IF=xy_IF,n_avg=avg_times)
        return False, {}


if __name__ == "__main__":
    
    """ Fill in """
    execution:bool = 1
    chip_info_restore:bool = 0
    DRandIP = {"dr":"drke","last_ip":"242"}
    ro_elements = ['q0']
    couplers = []

    
    """ Optional paras """

    freq_pts:int = 200
    freq_span_Hz:float = 400e6
    sweet_flux_shifter:float = 0
    xy_IF = 50e6
    avg_n:int = 300
    

    ro_p_max:float = 0.26#the output figure will show (ro_p_max)**2
    p_pts = 40




    """ Preparations """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    if ro_elements == 'all':#?
        ro_elements = list(Fctrl.keys())#?
    chip_info = cds.Chip_file(QD_agent=QD_agent)


    """ Running """
    FQ_results = {}
    check_again =[]
    Cctrl = coupler_zctrl(DRandIP["dr"],cluster,QD_agent.Fluxmanager.build_Cctrl_instructions(couplers,'i'))
    for qubit in ro_elements:
        if not QD_agent.Fluxmanager.get_offsweetspot_button(qubit):#?
        # if QD_agent.Fluxmanager.get_offsweetspot_button(qubit): #**when at off-sweet spot**
            init_system_atte(QD_agent.quantum_device,list([qubit]),ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'ro'),xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit,'xy'))
            # Cctrl['c0'](-0.15)
            # Cctrl['c1'](0.104)
            
            nc_path,trustable = StarkShift_executor(QD_agent, meas_ctrl,ro_elements, max_p=ro_p_max, run=execution, p_pts=p_pts, fpts=freq_pts, f_sapn_Hz=freq_span_Hz, avg_times=avg_n, xy_IF=xy_IF)
            #,new_ans
            
            # Cctrl['c0'](0)
            # Cctrl['c1'](0)
            cluster.reset()

            """ Storing """
            if  trustable:
                update_by_StarkShift(QD_agent,new_ans,qubit)
                QD_agent.QD_keeper()
                if chip_info_restore:
                    chip_info.update_by_StarkShift(qb=qubit, result=new_ans)
            else:
                check_again.append(qubit)    

    """ Close """
    print('Stark shift measurement done!')
    # if len(check_again) != 0:
    #     warning_print(f"qubits to check again: {check_again}")
    shut_down(cluster,Fctrl,Cctrl)

  
    


