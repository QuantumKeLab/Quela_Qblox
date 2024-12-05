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


def Cavity_powerDep_spec(quantum_device:QuantumDevice,ro_bare_guess:dict,ro_span_Hz:int=5e6,ro_p_min:float=0.1,ro_p_max:float=0.7,n_avg:int=1000,f_points:int=200,p_points:int=200,run:bool=True,q:str='q0',Save:bool=False):
    hw_c= quantum_device.hardware_config()
    sche_func = One_tone_sche
        
    analysis_result = {}
    analysis_result[q]= []
    data_save = {}
    ro_f_center = ro_bare_guess[q]
    ro_f_samples = linspace(ro_f_center-ro_span_Hz,ro_f_center+ro_span_Hz,f_points)
    ro_p_samples = linspace(ro_p_min,ro_p_max,p_points)
    freq = ManualParameter(name="freq", unit="Hz", label="Frequency")
    freq.batched = True
    
    ro_pulse_amp = ManualParameter(name="ro_amp", unit="", label="Readout pulse amplitude")
    ro_pulse_amp.batched = False
    
    
    spec_sched_kwargs = dict(   
        frequencies=freq,
        q=q,
        R_amp=ro_pulse_amp,
        R_duration=R_duration,
        R_integration=R_integration,
        R_inte_delay=R_inte_delay,
        powerDep=True,
    )
    exp_kwargs= dict(sweep_F=['start '+'%E' %ro_f_samples[0],'end '+'%E' %ro_f_samples[-1]],
                     f_points=f_points,
                     Power=['start '+'%E' %ro_p_samples[0],'end '+'%E' %ro_p_samples[-1]],
                     p_points=p_points,)
    
    if run:
        if Save is True:
            write_presave_txt(Save_filepath+q+'_'+'One-tone-powerDep_')
            
        else: pass
    
        gettable = ScheduleGettable(
            quantum_device,
            schedule_function=sche_func, 
            schedule_kwargs=spec_sched_kwargs,
            real_imag=False,
            batched=True,
        )
        quantum_device.cfg_sched_repetitions(n_avg)
        meas_ctrl.gettables(gettable)
        meas_ctrl.settables([freq,ro_pulse_amp])
        meas_ctrl.setpoints_grid((ro_f_samples,ro_p_samples))
        
        
        
        rp_ds = meas_ctrl.run("One-tone-powerDep")
        rp_ds
        Amp,phase= dataset_to_array(dataset=rp_ds,dims=2)
        data= [Amp.transpose(),phase.transpose()]
        for i in range(p_points):
            data_fit = Notch_type_resonator_fit(ro_f_samples,data[0][i],data[1][i],target_Ql=None)
            analysis_result[q].append(data_fit)
        data_save[q]= [ro_f_samples,ro_p_samples,data]
        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q,globals()))
        
        if Save is True:
            write_txt(Save_filepath+q+'_'+'One-tone-powerDep', "One_tone_powerDep_kwargs: Meas.qubit="+q,'\n\n',exp_kwargs,'\n\n',Experi_info(q,globals()),'\n\n','avg_times='+str(n_avg),'\n\n',hw_c)
            save_exper_info(Save_filepath+q+'_'+'One-tone-powerDep_'+"Exp_parameters_save_",globals())
            save_data(Save_filepath+q+'_'+'One-tone-powerDep',locals(),'rp_ds')
        else: pass
        
    else:
        sweep_para1= np.array(ro_f_samples[:n_s])
        sweep_para2= np.array(ro_p_samples[:2])
        spec_sched_kwargs['frequencies']= sweep_para1.reshape(sweep_para1.shape or (1,))
        spec_sched_kwargs['R_amp']= {q:sweep_para2.reshape(sweep_para2.shape or (1,))[0]}
        pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
        show_args(exp_kwargs, title="One_tone_powerDep_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q,globals()))
    return analysis_result,data_save


execute = True

if execute:
    PD_results,data_bank_pd = Cavity_powerDep_spec(quantum_device,
                               ro_bare_guess=ro_bare,
                               ro_span_Hz=8e6,
                               ro_p_min=0.05,
                               ro_p_max=1,
                               n_avg=100,
                               f_points=101,
                               p_points=11,
                               run=True,
                               q='q0',
                               Save=False,)



def One_tone_sche(
    frequencies: np.ndarray,
    q:str,
    R_amp: dict,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    powerDep:bool,
    repetitions:int=1,    
) -> Schedule:
    
    sched = Schedule("One tone spectroscopy (NCO sweep)",repetitions=repetitions)
    sched.add_resource(ClockResource(name=q+ ".ro", freq=frequencies.flat[0]))
    
    for acq_idx, freq in enumerate(frequencies):
        
        sched.add(SetClockFrequency(clock= q+ ".ro", clock_freq_new=freq))
        sched.add(Reset(q))
        sched.add(IdlePulse(duration=5000*1e-9), label=f"buffer {acq_idx}")
        spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=powerDep)
        
        Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_index=acq_idx,acq_channel=0,single_shot=False,get_trace=False,trace_recordlength=0)
     
    return sched



def Readout_amp_optimization(quantum_device:QuantumDevice,R_amp_low:float,R_amp_high:float,shots:int=1000,points:int=200,run:bool=True,q:str='q0',Save:bool=False):
    sche_func= Qubit_amp_SS_sche
    LO= f01[q]+IF
    hw_c=set_LO_frequency(quantum_device,q=q,module_type='drive',LO_frequency=LO)    
    analysis_result = {}
    raw_data={}
    raw_data[q]=[]
    analysis_result[q]=[]
    amp_samples = linspace(R_amp_low,R_amp_high,points)
    exp_kwargs= dict(sweep_amp=['start '+'%E' %amp_samples[0],'end '+'%E' %amp_samples[-1]],
                     shots=shots,points=points)
    def state_dep_sched(ini_state,amp,data):
        
        spec_sched_kwargs = dict(  
            q= q,
            ini_state= ini_state,
            pi_amp=pi_amp,
            pi_Du=pi_Du,
            R_amp=amp,
            R_duration=R_duration,
            R_integration=R_integration,
            R_inte_delay=R_inte_delay
        )
        
        if run:
            gettable = ScheduleGettable(
                quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=spec_sched_kwargs,
                real_imag=True,
                batched=True,
            )
            
            quantum_device.cfg_sched_repetitions(shots)
            rs_ds = gettable.get()
            rs_ds
            data[ini_state] = rs_ds
        else:
            pulse_preview(quantum_device,sche_func,spec_sched_kwargs)
    if run:
        if Save is True:
            write_presave_txt(Save_filepath+q+'_'+'Readout_amp_cal_')
            
        else: pass 
    tau= R_integration[q] 
    if run:
        for i in range(points):
            data = {}
            state_dep_sched('g',amp_samples[i],data)
            state_dep_sched('e',amp_samples[i],data)
            print(i)
            raw_data[q].append(data)
        for i in range(points):     
            print(i)
            analysis_result[q].append(Qubit_state_single_shot_fit_analysis(raw_data[q][i],T1=T1,tau=tau,f01=f01[q]))   
            
        show_args(exp_kwargs, title="Readout_amp_cal_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q,globals()))
    else:        
        data = {}
        state_dep_sched('g',amp_samples[0],data)
        state_dep_sched('e',amp_samples[0],data)
        show_args(exp_kwargs, title="Readout_amp_cal_kwargs: Meas.qubit="+q)
        show_args(Experi_info(q,globals()))
    if Save is True:
        write_txt(Save_filepath+q+'_'+'Readout_amp_cal', "Readout_amp_cal_kwargs: Meas.qubit="+q,'\n\n',exp_kwargs,'\n\n',Experi_info(q,globals()),'\n\n',hw_c)
        save_exper_info(Save_filepath+q+'_'+'Readout_amp_cal_'+"Exp_parameters_save_",globals())
        save_data(Save_filepath+q+'_'+'Readout_amp_cal',locals(),'raw_data')
    else: pass 
    return [analysis_result, amp_samples],raw_data



execute = True

if execute:
    Readout_amp_cal_results,Readout_amp_cal_raw_data = Readout_amp_optimization(quantum_device,
                                                                               R_amp_low=0.01,
                                                                               R_amp_high=0.6,
                                                                               shots=10000,
                                                                               points=41,
                                                                               run=True,
                                                                               q='q0',
                                                                               Save=True,)





def Qubit_amp_SS_sche(
    q:str,
    ini_state:str,
    pi_amp: dict,
    pi_Du: dict,
    R_amp: any,
    R_duration: dict,
    R_integration:dict,
    R_inte_delay:float,
    repetitions:int=1,
) -> Schedule:

    sched = Schedule("Single shot", repetitions=repetitions)

    sched.add(Reset(q))
    
    sched.add(IdlePulse(duration=5000*1e-9))
    
    spec_pulse = Readout(sched,q,R_amp,R_duration,powerDep=True)
    
    if ini_state=='e': 
        X_pi_p(sched,pi_amp,pi_Du,q,spec_pulse,freeDu=0)
        
    else: None
    Integration(sched,q,R_inte_delay,R_integration,spec_pulse,acq_index=0,acq_channel=0,single_shot=True,get_trace=False,trace_recordlength=0)

    return sched

