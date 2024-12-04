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
from Modularize.support import QDmanager, Data_manager,init_system_atte, init_meas, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import Qubit_amp_SS_sche, set_LO_frequency, pulse_preview, Qubit_state_single_shot_fit_analysis


try:
    from qcat.analysis.state_discrimination.discriminator import train_GMModel # type: ignore
    from qcat.visualization.readout_fidelity import plot_readout_fidelity
    from qblox_drive_AS.analysis.OneShotAna import a_OSdata_analPlot
    mode = "AS"
except:
    mode = "WeiEn"


def Qubit_state_single_shot(QD_agent: QDmanager, shots: int = 1000, run: bool = True, q: str = 'q1', IF: float = 250e6,
                            Experi_info: dict = {}, ro_amp_values: list = [0.1], T1: float = 15e-6, exp_idx: int = 0,
                            parent_datafolder: str = '', plot: bool = False):
    data = {}
    analysis_result = {}

    for ro_amp in ro_amp_values:
        qubit_info = QD_agent.quantum_device.get_element(q)
        qubit_info.measure.integration_time(2e-6)
        qubit_info.measure.pulse_duration(2e-6)
        qubit_info.measure.pulse_amp(ro_amp)  # Set the readout amplitude directly to the given value

        print("Integration time ", qubit_info.measure.integration_time() * 1e6, "µs")
        print("Reset time ", qubit_info.reset.duration() * 1e6, "µs")
        print("RO amplitude: ", qubit_info.measure.pulse_amp())

        sche_func = Qubit_amp_SS_sche
        LO = qubit_info.clock_freqs.f01() + IF
        set_LO_frequency(QD_agent.quantum_device, q=q, module_type='drive', LO_frequency=LO)

        exp_kwargs = dict(shots=shots)

        def state_dep_sched(ini_state: str):
            print(f"Shotting for |{ini_state}> with readout amplitude {ro_amp}")
            sched_kwargs = dict(
                q=q,
                ini_state=ini_state,
                pi_amp={str(q): qubit_info.rxy.amp180() * 1},
                pi_dura={str(q): qubit_info.rxy.duration()},
                R_amp={str(q): ro_amp},  # Use the current readout amplitude value
                R_duration={str(q): qubit_info.measure.pulse_duration()},
                R_integration={str(q): qubit_info.measure.integration_time()},
                R_inte_delay=qubit_info.measure.acq_delay(),
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
                data[(ini_state, ro_amp)] = ss_ds
                show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit=" + q)
                if Experi_info != {}:
                    show_args(Experi_info(q))

            else:
                pulse_preview(QD_agent.quantum_device, sche_func, sched_kwargs)
                show_args(exp_kwargs, title="Single_shot_kwargs: Meas.qubit=" + q)
                if Experi_info != {}:
                    show_args(Experi_info(q))

        tau = qubit_info.measure.integration_time()
        state_dep_sched('g')
        state_dep_sched('e')

    SS_dict = {f"g_{ro_amp}": {"dims": ("I", "Q"), "data": array(data[('g', ro_amp)])} for ro_amp in ro_amp_values}
    SS_dict.update({f"e_{ro_amp}": {"dims": ("I", "Q"), "data": array(data[('e', ro_amp)])} for ro_amp in ro_amp_values})

    SS_ds = Dataset.from_dict(SS_dict)
    nc_path = Data_manager().save_raw_data(QD_agent=QD_agent, ds=SS_ds, qb=q, exp_type='ss', label=exp_idx,
                                           specific_dataFolder=parent_datafolder, get_data_loc=True)

    if mode == "WeiEn" and plot:
        print("Under built-in analysis...")
        analysis_result[q] = Qubit_state_single_shot_fit_analysis(data, T1=T1, tau=tau)
    else:
        analysis_result[q] = []

    return analysis_result, nc_path


def SS_executor(QD_agent: QDmanager, cluster: Cluster, Fctrl: dict, target_q: str, shots: int = 10000,
                execution: bool = True, data_folder='', plot: bool = True, ro_amp_values: list = [0.1], 
                exp_label: int = 0, save_every_pic: bool = False, IF: float = 250e6):

    # Set proper Z-bias for the target qubit
    Fctrl[target_q](float(QD_agent.Fluxmanager.get_proper_zbiasFor(target_q)))

    # Perform Single Shot Measurement with different readout power values
    SS_result, nc = Qubit_state_single_shot(
        QD_agent,
        shots=shots,
        ro_amp_values=ro_amp_values,  # Pass list of readout power values
        run=execution,
        q=target_q,
        IF=IF,
        plot=plot,
        parent_datafolder=data_folder,
        exp_idx=exp_label
    )

    # Reset the flux for the target qubit
    Fctrl[target_q](0.0)
    cluster.reset()

    # Analysis part
    if mode == "WeiEn":
        if plot:
            Qubit_state_single_shot_plot(SS_result[target_q], Plot_type='both', y_scale='log')
            effT_mk, ro_fidelity, thermal_p = 0, 0, 0
        else:
            effT_mk, ro_fidelity, thermal_p = 0, 0, 0
    else:
        # Use the analysis function to extract necessary data
        thermal_p, effT_mk, ro_fidelity, p10_precentage, snr, power_snr_dB = a_OSdata_analPlot(
            QD_agent, target_q, nc, plot, save_pic=save_every_pic
        )

    return thermal_p, effT_mk, ro_fidelity

if __name__ == '__main__':

    """ Fill in """
    execute: bool = 1
    repeat: int = 1
    DRandIP = {"dr": "dr4", "last_ip": "81"}
    ro_elements = {'q1': {"roAmp_factor": 1}}
    couplers = []

    """ Optional paras (don't use is better) """
    ro_atte_degrade_dB: int = 0  # multiple of 2
    shot_num: int = 30000
    xy_IF = 250e6

    # Define multiple readout amplitudes (Voltages) to iterate over
    ro_amp_values = [0.01, 0.1, 0.3, 0.6]

    """ Iteration """
    snr_rec, effT_rec, thermal_pop = {}, {}, {}
    for qubit in ro_elements:
        for i in range(repeat):
            start_time = time.time()

            """ Preparation """
            slightly_print(f"The {i}th OS:")
            QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"], ip_label=DRandIP["last_ip"])
            QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path, mode='l')
            QD_agent.Notewriter.modify_DigiAtte_For(-ro_atte_degrade_dB, qubit, 'ro')

            """ Running """
            Cctrl = coupler_zctrl(DRandIP["dr"], cluster, QD_agent.Fluxmanager.build_Cctrl_instructions(couplers, 'i'))
            if i == 0:
                snr_rec[qubit], effT_rec[qubit], thermal_pop[qubit] = [], [], []
            init_system_atte(QD_agent.quantum_device, list([qubit]), xy_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit, 'xy'), ro_out_att=QD_agent.Notewriter.get_DigiAtteFor(qubit, 'ro'))
            
            # Use multiple readout amplitudes
            info = SS_executor(
                QD_agent, cluster, Fctrl, qubit, execution=execute, shots=shot_num,
                ro_amp_values=ro_amp_values, plot=True if repeat == 1 else False, exp_label=i, IF=xy_IF
            )

            snr_rec[qubit].append(info[2])
            effT_rec[qubit].append(info[1])
            thermal_pop[qubit].append(info[0] * 100)

            """ Storing """ 
            if execute and repeat == 1:
                keep = 'y' if ro_amp_scaling == 1 and ro_atte_degrade_dB == 0 else mark_input(f"Keep this RO amp for {qubit}?[y/n]")
                if keep.lower() in ['y', 'yes']:
                    QD_agent.QD_keeper()

            """ Close """    
            shut_down(cluster, Fctrl, Cctrl)
            end_time = time.time()
            slightly_print(f"time cost: {round(end_time - start_time, 1)} secs")

    for qubit in effT_rec:
        highlight_print(f"{qubit}: {round(median(array(effT_rec[qubit])), 2)} +/- {round(std(array(effT_rec[qubit])), 3)} mK")

        Data_manager().save_histo_pic(QD_agent, effT_rec, qubit, mode="ss")
        Data_manager().save_histo_pic(QD_agent, thermal_pop, qubit, mode="pop")

        
        

        
    
