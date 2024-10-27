import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity
from qcat.visualization.readout_fidelity import plot_readout_fidelity
from xarray import Dataset, open_dataset
from Modularize.analysis.Radiator.RadiatorSetAna import OSdata_arranger
from Modularize.support.UserFriend import *
from Modularize.support.QDmanager import QDmanager, Data_manager
from numpy import array, moveaxis, mean, std, median, arange
import matplotlib.pyplot as plt
from Modularize.analysis.Radiator.RadiatorSetAna import sort_files
from qcat.analysis.state_discrimination import p01_to_Teff
import numpy as np
from matplotlib.ticker import MultipleLocator
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set
from Modularize.analysis.OneShotAna import a_OSdata_analPlot
import re


if __name__=="__main__":

    """Fill in"""
    QD_agent_path=r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_10_25\DRKE#242_SumInfo.pkl"
    target_q='q0'
    data_folder = r'C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_25\15mK_IntegrationTime_1250ns\SS'
    show_each_plot=False

    """Multifiles Analysis"""
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()
    results = {'time': [], 'p01': [], 'effT_mK': []}

    # 遍历 data_folder 中的文件
    for file in os.listdir(data_folder):
        if file.endswith('.nc'):  # 假设文件具有 .nc 扩展名
            nc_path = os.path.join(data_folder, file)
            p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, target_q, nc_path,plot=show_each_plot)
            
            # 提取时间信息
            match = re.search(r'H(\d+)M(\d+)S(\d+)', file)
            if match:
                hours = int(match.group(1))
                minutes = int(match.group(2))
                seconds = int(match.group(3))
                total_time_seconds = hours * 3600 + minutes * 60 + seconds
                results['time'].append(total_time_seconds)
            
            # 检查 p01 和 effT_mK 是否为有效值
            if not (np.isnan(p01) or p01 < 0):
                results['p01'].append(p01)
            
            if not (np.isnan(effT_mK) or effT_mK < 0):
                results['effT_mK'].append(effT_mK)

# 创建 result_dict 以便用于保存图形
result_dict = {q: np.array(v) for q, v in results.items()}

# 绘制时间对 p01 和 effT_mK 的图
time_array = result_dict['time']
p01_array = result_dict['p01']
effT_mK_array = result_dict['effT_mK']

# 确保数组的长度一致
min_length = min(len(time_array), len(p01_array), len(effT_mK_array))
time_array = time_array[:min_length]
p01_array = p01_array[:min_length]
effT_mK_array = effT_mK_array[:min_length]

# 绘制图形
plt.figure(figsize=(10, 6))
plt.plot(time_array, p01_array, label='p01', marker='o', linestyle='-', color='blue')
plt.plot(time_array, effT_mK_array, label='effT_mK', marker='x', linestyle='-', color='orange')

plt.title('Time vs. p01 and effT_mK')
plt.xlabel('Time (seconds)')
plt.ylabel('Values')
plt.legend()
plt.grid()
plt.tight_layout()

# 显示图形
plt.show()

# 保存 histogram 图
Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder='path_to_save_thermalpop_plots')
Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder='path_to_save_effT_mK_plots')
