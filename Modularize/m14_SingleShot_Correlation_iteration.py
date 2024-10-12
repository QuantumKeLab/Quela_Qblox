import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from Modularize.support.UserFriend import *
from qcodes.parameters import ManualParameter
from numpy import array, linspace, median, std
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support.Pulse_schedule_library import Qubit_state_single_shot_plot
from Modularize.support import QDmanager, Data_manager,init_system_atte, init_meas, shut_down, coupler_zctrl
from Modularize.support.Pulse_schedule_library import Qubit_SS_sche, Qubit_SS_Correlation_sche, set_LO_frequency, pulse_preview, Qubit_state_single_shot_fit_analysis

try:
    from qcat.analysis.state_discrimination.discriminator import train_GMModel # type: ignore
    from qcat.visualization.readout_fidelity import plot_readout_fidelity
    from Modularize.analysis.OneShotAna import a_OSdata_analPlot, a_OSdata_correlation_analPlot
    mode = "AS"
except:
    mode = "WeiEn"

def Qubit_state_single_shot(QD_agent:QDmanager, shots:int=1000, run:bool=True, q:str='q1', IF:float=250e6, Experi_info:dict={}, ro_amp_factor:float=1, T1:float=15e-6, exp_idx:int=0, parent_datafolder:str='', plot:bool=False, correlate_delay:float=0.0):
    qubit_info = QD_agent.quantum_device.get_element(q)
    qubit_info.measure.integration_time(1.2e-6)
    qubit_info.measure.pulse_duration(1.2e-6)
    print("Integration time ",qubit_info.measure.integration_time()*1e6, "µs")
    print("Reset time ", qubit_info.reset.duration()*1e6, "µs")    
    sche_func=Qubit_SS_Correlation_sche
    LO= qubit_info.clock_freqs.f01() + IF
    set_LO_frequency(QD_agent.quantum_device, q=q, module_type='drive', LO_frequency=LO)
    if ro_amp_factor != 1:
        qubit_info.measure.pulse_amp(ro_amp_factor*qubit_info.measure.pulse_amp())
        eyeson_print(f"The new RO amp = {round(qubit_info.measure.pulse_amp(),2)}")
    else:
        eyeson_print(f"RO amp = {qubit_info.measure.pulse_amp()}")
    data = {}
    analysis_result = {}
    exp_kwargs = dict(shots=shots)
    print(qubit_info.rxy.amp180())  
    
    def state_dep_sched(ini_state:str):
        slightly_print(f"Shotting for |{ini_state}> with correlate_delay={correlate_delay} µs")
        sched_kwargs = dict(   
            q=q,
            ini_state=ini_state,
            pi_amp={str(q):qubit_info.rxy.amp180()},
            pi_dura={str(q):qubit_info.rxy.duration()},
            R_amp={str(q):qubit_info.measure.pulse_amp()},
            R_duration={str(q):qubit_info.measure.pulse_duration()},
            R_integration={str(q):qubit_info.measure.integration_time()},
            R_inte_delay=qubit_info.measure.acq_delay(),
            correlate_delay=correlate_delay * 1e-6  # Set correlate_delay in seconds
        )

        if run:
            gettable = ScheduleGettable(
                QD_agent.quantum_device,
                schedule_function=sche_func, 
                schedule_kwargs=sched_kwargs,
                real_imag=True,
                batched=True,
            )
            QD_agent.quantum_device.cfg_sched_repetitions(shots)
            ss_ds = gettable.get()
            print(array(ss_ds).shape)
            data[ini_state] = ss_ds
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit=" + q)
            if Experi_info != {}:
                show_args(Experi_info(q))
        else:
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)
            
            show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit="+q)
            if Experi_info != {}:
                show_args(Experi_info(q))

    tau = qubit_info.measure.integration_time()        
    state_dep_sched('g')

    SS_dict = {
        "g": {"dims": ("I", "Q"), "data": array(data['g'])},
    }

    SS_ds = Dataset.from_dict(SS_dict)
    nc_path = Data_manager().save_raw_data(QD_agent=QD_agent, ds=SS_ds, qb=q, exp_type='ss', label=exp_idx, specific_dataFolder=parent_datafolder, get_data_loc=True)
    
    if mode == "WeiEn" and plot: 
        slightly_print("Under built-in analysis...")
        analysis_result[q] = Qubit_state_single_shot_fit_analysis(data, T1=T1, tau=tau) 
    else:
        analysis_result[q] = []

    return analysis_result, nc_path


def SS_executor(QD_agent:QDmanager, cluster:Cluster, Fctrl:dict, target_q:str, shots:int=10000, execution:bool=True, data_folder='', roAmp_modifier:float=1, IF:float=250e6):
    Fctrl[target_q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(target_q)))

    for i, correlate_delay in enumerate(range(0,31)):  # range(0, 11) is 0 to 10 µs
        plot = (i == 0)  # Only plot the first measurement with correlate_delay=0
        nc = Qubit_state_single_shot(QD_agent,
                                    shots=shots,
                                    run=execution,
                                    q=target_q,
                                    parent_datafolder=data_folder,
                                    ro_amp_factor=roAmp_modifier,
                                    exp_idx=i,
                                    plot=plot,
                                    IF=IF,
                                    correlate_delay=correlate_delay)
        
    Fctrl[target_q](0.0)
    cluster.reset()

    if plot:
        a_OSdata_correlation_analPlot(nc,plot)

if __name__ == '__main__':
    execute = True
    DRandIP = {"dr": "drke", "last_ip": "242"}
    ro_elements = {'q1': {"roAmp_factor": 1}}
    couplers = []


    """ Optional paras (don't use is better) """
    ro_atte_degrade_dB:int = 0 # multiple of 2 
    shot_num:int = 50000
    xy_IF = 250e6


    for qubit in ro_elements:
        start_time = time.time()
        QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"], ip_label=DRandIP["last_ip"])
        QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path, mode='l')
        QD_agent.Notewriter.modify_DigiAtte_For(-ro_atte_degrade_dB, qubit, 'ro')

        Cctrl = coupler_zctrl(DRandIP["dr"], cluster, QD_agent.Fluxmanager.build_Cctrl_instructions(couplers, 'i'))
        init_system_atte(QD_agent.quantum_device, [qubit], 
                         xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit, 'xy'),
                         ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit, 'ro'))
        ro_amp_scaling = ro_elements[qubit]["roAmp_factor"]
        SS_executor(QD_agent, cluster, Fctrl, qubit, execution=execute, shots=shot_num, roAmp_modifier=ro_amp_scaling, IF=xy_IF)
        
        if ro_amp_scaling !=1 or ro_atte_degrade_dB != 0:
            keep = mark_input(f"Keep this RO amp for {qubit}?[y/n]")
        else:
            keep = 'y'

        """ Storing """ 
        if execute:
            if keep.lower() in ['y', 'yes']:
                QD_agent.QD_keeper() 

        """ Close """    
        shut_down(cluster,Fctrl,Cctrl)
        end_time = time.time()
        slightly_print(f"time cose: {round(end_time-start_time,1)} secs")
