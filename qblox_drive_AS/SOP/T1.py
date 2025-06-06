import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from xarray import Dataset
from numpy import array, arange
from qblox_drive_AS.support.UserFriend import *
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support import Data_manager
from quantify_scheduler.gettables import ScheduleGettable
from qblox_drive_AS.support import compose_para_for_multiplexing
from qblox_drive_AS.support.Pulse_schedule_library import pulse_preview
from qblox_drive_AS.support.Pulser import ScheduleConductor
from qblox_drive_AS.support.Pulse_schedule_library import Schedule, Readout, Multi_Readout, Integration, electrical_delay
from quantify_scheduler.operations.gate_library import Reset


class EnergyRelaxPS(ScheduleConductor):
    def __init__(self):
        super().__init__()
        self._time_samples:dict = {}
        self._os_mode:bool = False
        self._repeat:int = 1
        self._avg_n:int = 300

    @property
    def time_samples( self ):
        return self._time_samples
    @time_samples.setter
    def set_time_samples(self, time_samples:dict):
        if not isinstance(time_samples,dict):
            raise TypeError("Arg 'time_samples' must be a dict !")
        self._time_samples = time_samples
    @property
    def repeat( self ):
        return self._repeat
    @repeat.setter
    def set_repeat(self, repeat_num:int):
        if not isinstance(repeat_num,(int,float)):
            raise TypeError("Ard 'repeat_num' must be a int or float !")
        self._repeat = int(repeat_num)
    @property
    def os_mode( self ):
        return self._os_mode
    @os_mode.setter
    def set_os_mode(self, os_mode:bool):
        if not isinstance(os_mode,bool):
            if os_mode in [0,1]:
                pass
            else:
                raise TypeError("Arg 'os_mode' must be a bool, 0 or 1 !")
        
        self._os_mode = os_mode
    @property
    def n_avg( self ):
        return self._avg_n
    @n_avg.setter
    def set_n_avg(self, avg_num:int):
        if not isinstance(avg_num,(int,float)):
            raise TypeError("Ard 'avg_num' must be a int or float !")
        self._avg_n = int(avg_num)


    def __PulseSchedule__(self, 
        freeduration:dict,
        pi_amp: dict,
        pi_dura:dict,
        R_amp: dict,
        R_duration: dict,
        R_integration:dict,
        R_inte_delay:dict,
        repetitions:int=1,
        singleshot:bool=False 
        )->Schedule:

        qubits2read = list(freeduration.keys())
        sameple_idx = array(freeduration[qubits2read[0]]).shape[0]
        sched = Schedule("T1", repetitions=repetitions)
        for acq_idx in range(sameple_idx):
            for qubit_idx, q in enumerate(qubits2read):
                freeDu = freeduration[q][acq_idx]
                sched.add(Reset(q))
            
                if qubit_idx == 0:
                    spec_pulse = Readout(sched,q,R_amp,R_duration)
                else:
                    Multi_Readout(sched,q,spec_pulse,R_amp,R_duration)
            
            
                self.QD_agent.Waveformer.X_pi_p(sched,pi_amp,q,pi_dura[q],spec_pulse,freeDu+electrical_delay)
                
                Integration(sched,q,R_inte_delay[q],R_integration,spec_pulse,acq_idx,acq_channel=qubit_idx,single_shot=singleshot,get_trace=False,trace_recordlength=0)
        
        self.schedule = sched
        return sched

    def __SetParameters__(self, *args, **kwargs):
        
        self.__time_data_idx = arange(len(list(self._time_samples.values())[0]))
        for q in self._time_samples:
            qubit_info = self.QD_agent.quantum_device.get_element(q)
            if max(self._time_samples[q]) >= qubit_info.reset.duration():
                qubit_info.reset.duration(qubit_info.reset.duration()+50e-6)
                eyeson_print(f"{q} Reset time had been extended to: {round(qubit_info.reset.duration()*1e6,0)} µs")
            else:
                eyeson_print(f"{q} Reset time: {round(qubit_info.reset.duration()*1e6,0)} µs")
        
        self.__Para_free_Du = ManualParameter(name="free_Duration", unit="s", label="Time")
        self.__Para_free_Du.batched = True
        self.__Para_repeat = ManualParameter(name="repeat", unit="n", label="Count")
        self.__Para_repeat.batched = False
        self.__repeat_data_idx = arange(self._repeat)

        if self._os_mode:
            self.__one_shot_para =  ManualParameter(name="Shot")

        self.__sched_kwargs = dict(
        freeduration=self._time_samples,
        singleshot=self._os_mode,
        pi_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d1'),
        pi_dura=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'d3'),
        R_amp=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r1'),
        R_duration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r3'),
        R_integration=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r4'),
        R_inte_delay=compose_para_for_multiplexing(self.QD_agent,self._time_samples,'r2'),
        )
    
    def __Compose__(self, *args, **kwargs):
        
        if self._execution:
            self.__gettable = ScheduleGettable(
                self.QD_agent.quantum_device,
                schedule_function=self.__PulseSchedule__,
                schedule_kwargs=self.__sched_kwargs,
                real_imag=True,
                batched=True,
                num_channels=len(list(self._time_samples.keys())),
                )
            self.QD_agent.quantum_device.cfg_sched_repetitions(self._avg_n)
            self.meas_ctrl.gettables(self.__gettable)
            if not self._os_mode:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,self.__repeat_data_idx))
            else:
                self.meas_ctrl.settables([self.__Para_free_Du,self.__one_shot_para,self.__Para_repeat])
                self.meas_ctrl.setpoints_grid((self.__time_data_idx,arange(self._avg_n),self.__repeat_data_idx))
        else:
            preview_para = {}
            for q in self._time_samples:
                preview_para[q] = array([self._time_samples[q][0],self._time_samples[q][-1]])
            self.__sched_kwargs['freeduration']= preview_para
    
    def __RunAndGet__(self, *args, **kwargs):
        
        if self._execution:
            ds = self.meas_ctrl.run('T1')
            dict_ = {}
            if not self._os_mode:
                for q_idx, q in enumerate(self._time_samples):
                    i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._time_samples[q].shape[0])
                    q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._time_samples[q].shape[0])
                    dict_[q] = (["mixer","repeat","idx"],array([i_data,q_data]))
                    time_values = list(self._time_samples[q])*2*self._repeat
                    dict_[f"{q}_x"] = (["mixer","repeat","idx"],array(time_values).reshape(2,self._repeat,self._time_samples[q].shape[0]))
                
                dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":self.__repeat_data_idx,"idx":self.__time_data_idx})
            else:
                dict_ = {}
                for q_idx, q in enumerate(self._time_samples):
                    i_data = array(ds[f'y{2*q_idx}']).reshape(self._repeat,self._avg_n,self._time_samples[q].shape[0])
                    q_data = array(ds[f'y{2*q_idx+1}']).reshape(self._repeat,self._avg_n,self._time_samples[q].shape[0])
                    dict_[q] = (["mixer","prepared_state","repeat","index","time_idx"],array([[i_data],[q_data]]))
                    time_values = list(self._time_samples[q])*2*self._repeat*self._avg_n
                    dict_[f"{q}_x"] = (["mixer","prepared_state","repeat","index","time_idx"],array(time_values).reshape(2,1,self._repeat,self._avg_n,self._time_samples[q].shape[0]))

                dataset = Dataset(dict_,coords={"mixer":array(["I","Q"]),"repeat":self.__repeat_data_idx,"prepared_state":array([1]),"index":arange(self._avg_n),"time_idx":self.__time_data_idx})
                
                
            dataset.attrs["execution_time"] = Data_manager().get_time_now()
            dataset.attrs["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            dataset.attrs["method"] = "Shot" if self._os_mode else "Average"
            dataset.attrs["system"] = "qblox"
            self.dataset = dataset

        else:
            pulse_preview(self.QD_agent.quantum_device,self.__PulseSchedule__,self.__sched_kwargs)


if __name__ == "__main__":
    t1 = EnergyRelaxPS()
    ps = t1.get_adjsutable_paras(display=True)