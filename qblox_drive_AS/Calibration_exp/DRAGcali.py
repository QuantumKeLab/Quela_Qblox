import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from xarray import Dataset
from numpy import array, arange, moveaxis
from qblox_drive_AS.support import QDmanager, Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import drag_coef_cali, pulse_preview


def drag_cali(QD_agent:QDmanager,meas_ctrl:MeasurementControl, drag_samples:dict, n_avg:int=300,run:bool=True):
    results = {}
    
    sche_func= drag_coef_cali
    
    dataset_2_nc = ""
    for q in drag_samples:
        results[q], results[f"{q}_dragcoef"] = [], []
        data_sample_idx = arange(drag_samples[q].shape[0])
        qubit_info = QD_agent.quantum_device.get_element(q)
        eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")

    
    Sweep_para = ManualParameter(name="drag_coef")
    Sweep_para.batched = True
    
    def operation_dep_exe(operation_idx:int)->dict:
        dict_ = {}
        sched_kwargs = dict(
            seq_combin_idx=operation_idx,
            drag_ratios=drag_samples,
            waveformer=QD_agent.Waveformer,
            pi_amp=compose_para_for_multiplexing(QD_agent,drag_samples,'d1'),
            XY_duration=compose_para_for_multiplexing(QD_agent,drag_samples,'d3'),
            R_amp=compose_para_for_multiplexing(QD_agent,drag_samples,'r1'),
            R_duration=compose_para_for_multiplexing(QD_agent,drag_samples,'r3'),
            R_integration=compose_para_for_multiplexing(QD_agent,drag_samples,'r4'),
            R_inte_delay=compose_para_for_multiplexing(QD_agent,drag_samples,'r2'),
            )
        
        
        if run:
            gettable = ScheduleGettable(
            QD_agent.quantum_device,
            schedule_function=sche_func,
            schedule_kwargs=sched_kwargs,
            real_imag=True,
            batched=True,
            num_channels=len(list(drag_samples.keys())),
            )
            
    
            QD_agent.quantum_device.cfg_sched_repetitions(n_avg)
            meas_ctrl.gettables(gettable)
            meas_ctrl.settables(Sweep_para)
            meas_ctrl.setpoints(data_sample_idx)
        
        
            ds = meas_ctrl.run("drag coef calibration")
            
            for q_idx, q in enumerate(drag_samples):
                I_data, Q_data = array(ds[f"y{2*q_idx}"]).tolist(), array(ds[f"y{2*q_idx+1}"]).tolist()
                dict_[q] = [I_data,Q_data] # shape in (mixer, dragcoef)
                dict_[f"{q}_dragcoef"] = [list(drag_samples[q])]*2
        else:
            preview_para = {}
            for q in drag_samples:
                preview_para[q] = array([drag_samples[q][0],drag_samples[q][-1]])
            sched_kwargs['drag_ratios']= preview_para
            pulse_preview(QD_agent.quantum_device,sche_func,sched_kwargs)

        return dict_
    
    for operation_idx, operations in enumerate(["(X,Y/2)","(Y,X/2)"]):
        slightly_print(f"operations: {operations} ...")
        dataDict = operation_dep_exe(operation_idx) # {"q0":[],"q0_dragcoef":[], ...}
        for var in dataDict:
            results[var].append(dataDict[var])
            if operation_idx == 1: results[var] = (["mixer", "operations", "dragCoef"],moveaxis(array(results[var]),0,1)) # shape (operations, mixer, dragCoef) -> (mixer, operations, dragCoef)
            
    
    if run:
        dataset_2_nc = Dataset(results,coords={"mixer":array(["I","Q"]),"operations":array(["(X,Y/2)","(Y,X/2)"]),"dragCoef":data_sample_idx})
        dataset_2_nc.attrs["execution_time"] = Data_manager().get_time_now()
        dataset_2_nc.attrs["method"] = "Average"
        dataset_2_nc.attrs["system"] = "qblox"
   
    return dataset_2_nc