import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import Dataset, open_dataset
from Modularize.support import QDmanager,Data_manager
import matplotlib.pyplot as plt
from numpy import array, median, std 
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis
from datetime import datetime 
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set
def time_label_sort(nc_file_name:str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0],"H%HM%MS%S")

folder = r"C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\T1_timeDep_q1_20241012"
QD_file_path = r"Modularize\QD_backup\2024_10_12\DRKE#242_SumInfo.pkl"
qs = ['q1']
time_cost = 30
# ref_t1 = 29.8
# ref_t1_err = 3.3

# import os, datetime
# def get_time_now():
#     """
#     Since we save the Xarray into netCDF, we use the current time to encode the file name.\n
#     Ex: 19:23:34 return H19M23S34 
#     """
#     current_time = datetime.datetime.now()
#     return f"H{current_time.hour}M{current_time.minute}S{current_time.second}"
'''Time label sort'''
## files = sorted([name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")],key=lambda name:time_label_sort(name))
'''Exp idx sort'''
files = [name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder,name)) and name.split(".")[-1] == "nc")]
sort_set(files,3)

QD_agent = QDmanager(QD_file_path)
QD_agent.QD_loader()
time = []
t1 = []
raw_data = []

for idx, file in enumerate(files) :
    path = os.path.join(folder,file)
    nc = open_dataset(path)
    samples = array(nc.x0)
    I,Q = dataset_to_array(dataset=nc,dims=1)
    data = IQ_data_dis(I,Q,ref_I=QD_agent.refIQ[qs[0]][0],ref_Q=QD_agent.refIQ[qs[0]][-1])*1000 # unit in mV
    raw_data.append(data)
   
    data_fit = T1_fit_analysis(data=data,freeDu=samples,T1_guess=25e-6)
    t1.append(data_fit.attrs['T1_fit']*1e6)
    time.append(time_cost*(idx+1))

Data_manager().save_histo_pic(QD_agent,{qs[0]:t1},qs[0],mode="t1")
raw_data = array(raw_data)
ref_t1 = median(array(t1))
ref_t1_err = std(array(t1))
plt.pcolormesh(array(time)/60,samples*1e6,raw_data.transpose(),cmap='RdBu')
plt.colorbar()
plt.plot(array(time)/60,t1,c='green')
plt.hlines(ref_t1+ref_t1_err,min(array(time))/60,max(array(time))/60,linestyles='--',colors='orange')
plt.hlines(ref_t1,min(array(time))/60,max(array(time))/60,linestyles='-',colors='orange')
plt.hlines(ref_t1-ref_t1_err,min(array(time))/60,max(array(time))/60,linestyles='--',colors='orange')
plt.xlabel("time past (min)")
plt.ylabel("Free evolution time (us)")
plt.title("Time dependent relaxation time")
plt.grid()
plt.tight_layout()
plt.show()