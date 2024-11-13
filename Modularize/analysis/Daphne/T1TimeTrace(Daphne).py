import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from xarray import Dataset, open_dataset
from Modularize.support import QDmanager, Data_manager
import matplotlib.pyplot as plt
from numpy import array, median, std
from Modularize.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T1_fit_analysis, Fit_analysis_plot
from datetime import datetime
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set
import numpy as np
 

def time_label_sort(nc_file_name: str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0], "H%HM%MS%S")

folder = r"C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_11_2\T1_timeDep"
QD_file_path = r"C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_11_2\DRKE#242_SumInfo.pkl"
T1_guess=100e-6
qs = ['q0']
time_cost = 50
time_mode = 'relative'  # 'relative' 或 'real'(not done yet)

# 確保目錄存在
"""P6F-2"""
output_dir = os.path.join("C:\\Users\\Ke Lab\\Documents\\GitHub\\Quela_Qblox\\Modularize", "Meas_raw", datetime.now().strftime("%Y_%m_%d"))
os.makedirs(output_dir, exist_ok=True)

"""P6F-1"""
# output_dir = os.path.join(r"C:\\Users\\admin\\Documents\\GitHub\\Quela_Qblox\\Quela_Qblox\\Modularize", "Meas_raw", datetime.now().strftime("%Y_%m_%d"))
# os.makedirs(output_dir, exist_ok=True)


"""Running"""
files = [name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder, name)) and name.split(".")[-1] == "nc")]
sort_set(files, 3)

QD_agent = QDmanager(QD_file_path)
QD_agent.QD_loader()
time = []
t1 = []
raw_data = []

for idx, file in enumerate(files):
    path = os.path.join(folder, file)
    nc = open_dataset(path)
    samples = array(nc.x0)
    I, Q = dataset_to_array(dataset=nc, dims=1)
    data = IQ_data_dis(I, Q, ref_I=QD_agent.refIQ[qs[0]][0], ref_Q=QD_agent.refIQ[qs[0]][-1]) * 1000  # unit in mV
    raw_data.append(data)
   
    try:
        data_fit = T1_fit_analysis(data=data, freeDu=samples, T1_guess=T1_guess)
        t1_fit = data_fit.attrs['T1_fit']
        if not np.isnan(t1_fit) and 0 < t1_fit < 50:
            t1.append(t1_fit * 1e6)
        else:
            if t1:  # Only replace if t1 already has valid entries
                mean_t1 = np.mean(t1)
                std_t1 = np.std(t1)
                t1.append(mean_t1 + std_t1)
    except ValueError as e:
        print(f"Error fitting data for file {file}: {e}")
        if t1:  # Only replace if t1 already has valid entries
            mean_t1 = np.mean(t1)
            std_t1 = np.std(t1)
            t1.append(mean_t1 + std_t1)
    
    time.append(time_cost * (idx + 1))


# 獲取實際時間
if time_mode == 'real':
    time_array = [datetime.strptime(file.split("_")[-1].split(".")[0], "H%HM%MS%S") for file in files]
    time_array = [t.strftime("%H:%M:%S") for t in time_array]  # 轉換為字符串格式
else:
    time_array = array(time) / 60

Data_manager().save_histo_pic(QD_agent, {qs[0]: t1}, qs[0], mode="t1")

raw_data = array(raw_data)

t1_hero = max(t1)
t1_hero_idx = t1.index(t1_hero)  # 找到最大 t1 值的索引
t1_hero_file = files[t1_hero_idx]  # 根據索引找到對應的 .nc 檔案名稱

print(t1)
print("Maximum t1 is:", t1_hero)
print("Corresponding .nc file:", t1_hero_file)

ref_t1 = median(array(t1))
ref_t1_err = std(array(t1))

plt.pcolormesh(time_array, samples * 1e6, raw_data.transpose(), cmap='RdBu')
plt.colorbar()
plt.plot(time_array, t1, c='green')
plt.hlines(ref_t1 + ref_t1_err, min(time_array), max(time_array), linestyles='--', colors='orange')
plt.hlines(ref_t1, min(time_array), max(time_array), linestyles='-', colors='orange')
plt.hlines(ref_t1 - ref_t1_err, min(time_array), max(time_array), linestyles='--', colors='orange')
plt.xlabel("Time past (min)" if time_mode == 'relative' else "Real time")
plt.ylabel("Free evolution time (us)")
plt.title("Time dependent relaxation time")
plt.grid()
plt.tight_layout()
plt.show()

