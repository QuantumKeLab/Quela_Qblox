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
from Modularize.support.Pulse_schedule_library import hist_plot

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
    
    for file in os.listdir(data_folder):
        if file.endswith('.nc'):  # Assuming the files have a .nc extension
            nc_path = os.path.join(data_folder, file)
            p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, target_q, nc_path,plot=show_each_plot)
            # 檢查 p01 和 effT_mK 是否為有效值
            if not (np.isnan(p01) or p01 < 0):  
                results['p01'].append(p01)
            if not (np.isnan(effT_mK) or effT_mK < 0): 
                results['effT_mK'].append(effT_mK)   

    print(results)
    # Convert results to a dictionary format suitable for plotting
    result_dict = {q: np.array(v) for q, v in results.items()}

Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=data_folder)
Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=data_folder)
print("Analysis done! Check out the histogram in folder.")



# # 假设这里定义你的hist_plot和QDmanager类等其他必要的函数和类

# def process_folders(base_folder):
#     results = {}
#     folder_names = ['40mK', '45mK', '50mK']
    
#     for folder in folder_names:
#         # 获取每个子文件夹的路径
#         ss_folder_path = os.path.join(base_folder, folder, 'SS')
        
#         if os.path.exists(ss_folder_path):
#             # 记录每个文件的 effT 值
#             effT_values = []
#             for file in os.listdir(ss_folder_path):
#                 if file.endswith('.nc'):  # 假设文件具有.nc扩展名
#                     nc_path = os.path.join(ss_folder_path, file)
#                     p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, 'q0', nc_path)
                    
#                     # 检查有效值
#                     if not (np.isnan(effT_mK) or effT_mK < 0):
#                         effT_values.append(effT_mK)

#                     # 绘制直方图
#                     hist_plot('effT_mK', {'effT_mK': effT_values}, title=folder.replace('mK', ''), save_path=os.path.join('path_to_save_histograms', f"{folder}_histogram.png"))

#             # 计算平均值和标准差
#             if effT_values:
#                 mean_effT = np.mean(effT_values)
#                 std_effT = np.std(effT_values)
#                 results[folder.replace('mK', '')] = (mean_effT, std_effT)
    
#     return results

# def plot_summary(results):
#     x_labels = list(results.keys())
#     means = [result[0] for result in results.values()]
#     std_devs = [result[1] for result in results.values()]
    
#     # 创建条形图
#     plt.figure(figsize=(8, 5))
#     plt.bar(x_labels, means, yerr=std_devs, capsize=5)
#     plt.xlabel('Temperature (mK)')
#     plt.ylabel('Average effT (mK)')
#     plt.title('Average effT with Standard Deviation')
#     plt.tight_layout()
#     plt.savefig('path_to_save_summary/summary_plot.png')
#     plt.show()

# # 使用示例
# base_folder = r'C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\2024_10_27'  # 资料夹A的路径
# results = process_folders(base_folder)
# plot_summary(results)



