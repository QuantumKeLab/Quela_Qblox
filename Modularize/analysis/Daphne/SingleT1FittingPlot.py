import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from xarray import Dataset, open_dataset
from Modularize.support import QDmanager
from numpy import array
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, Fit_analysis_plot


folder = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241031_DR5_FQV1+NCU-1\Qblox\FQV1_Q0_20241103\Q0\T1_timeDep"
QD_file_path = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241031_DR5_FQV1+NCU-1\Qblox\QD_backup\2024_11_3\DRKE#242_SumInfo.pkl"
T1_guess=100e-6
qs = ['q0']
file= r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241031_DR5_FQV1+NCU-1\Qblox\FQV1_Q0_20241103\Q0\T1_timeDep\DRKEq0_T1(557)_H7M35S47.nc"


QD_agent = QDmanager(QD_file_path)
QD_agent.QD_loader()


path = os.path.join(folder, file)
nc = open_dataset(path)
samples = array(nc.x0)
I, Q = dataset_to_array(dataset=nc, dims=1)
data = IQ_data_dis(I, Q, ref_I=QD_agent.refIQ[qs[0]][0], ref_Q=QD_agent.refIQ[qs[0]][-1]) * 1000  # unit in mV
   
   

data_fit = T1_fit_analysis(data=data, freeDu=samples, T1_guess=T1_guess)
Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)



 