import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from qblox_instruments import Cluster
from utils.tutorial_utils import show_args
from qcodes.parameters import ManualParameter
from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support import QDmanager, Data_manager, cds
from quantify_scheduler.gettables import ScheduleGettable
from numpy import std, arange, array, average, mean, ndarray
from quantify_core.measurement.control import MeasurementControl
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import init_meas, init_system_atte, shut_down, coupler_zctrl
from qblox_drive_AS.support.Pulse_schedule_library import Ramsey_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T2_fit_analysis, Fit_analysis_plot, Fit_T2_cali_analysis_plot, T1_fit_analysis
from xarray import open_dataset
path = r"Modularize\Meas_raw\2024_8_18\DRKEq4_T2(1)_H15M55S11.nc"
path1 = r"Modularize\Meas_raw\2024_8_18\DRKEq4_T2(1)_H16M10S39.nc"
ds = open_dataset(path)
time_samples = array(ds.x0)
I,Q= dataset_to_array(dataset=ds,dims=1)
data= IQ_data_dis(I,Q,ref_I=0,ref_Q=0)

data_fit = T2_fit_analysis(data=data,freeDu=time_samples,T2_guess=15e-6)
Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)

ds = open_dataset(path1)
time_samples = array(ds.x0)
I,Q= dataset_to_array(dataset=ds,dims=1)
data= IQ_data_dis(I,Q,ref_I=0,ref_Q=0)

data_fit = T2_fit_analysis(data=data,freeDu=time_samples,T2_guess=15e-6)
Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)
# T2_us[q] = data_fit.attrs['T2_fit']*1e6
# Real_detune[q] = data_fit.attrs['f']