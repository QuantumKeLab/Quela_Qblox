import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
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

"""One Folder, Multiple Files"""
# if __name__=="__main__":

#     """Fill in"""
#     QD_agent_path=r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_25\DRKE#242_SumInfo.pkl"
#     target_q='q0'
#     data_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\2024_10_25\15mK_IntegrationTime_1000ns\SS_30aveg"
#     show_each_plot=False

#     """Multifiles Analysis"""
#     Qmanager = QDmanager(QD_agent_path)
#     Qmanager.QD_loader()
#     results = {'p01': [], 'effT_mK': []}#'time': [], 
    
#     for file in os.listdir(data_folder):
#         if file.endswith('.nc'):  # Assuming the files have a .nc extension
#             nc_path = os.path.join(data_folder, file)
#             p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, target_q, nc_path,plot=show_each_plot)
#             # 檢查 p01 和 effT_mK 是否為有效值
#             if not (np.isnan(p01) or p01 < 0):  
#                 results['p01'].append(p01)
#             if not (np.isnan(effT_mK) or effT_mK < 0): 
#                 results['effT_mK'].append(effT_mK)   

#     print(results)
#     # Convert results to a dictionary format suitable for plotting
#     result_dict = {q: np.array(v) for q, v in results.items()}
    
# Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=data_folder)
# Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=data_folder)
# print("Analysis done! Check out the histogram in folder.")

"""Multiple Folders"""
def analyze_and_plot(Qmanager, target_q, folder_path, show_each_plot, save_folder):
    """針對特定資料夾中的 .nc 檔案進行分析並生成 histogram 圖"""
    results = {'p01': [], 'effT_mK': []}
    
    for file in os.listdir(folder_path):
        if file.endswith('.nc'):  # 確認檔案副檔名為 .nc
            nc_path = os.path.join(folder_path, file)
            p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, target_q, nc_path, plot=show_each_plot)
            
            # 檢查 p01 和 effT_mK 是否為有效值
            if not (np.isnan(p01) or p01 < 0):  
                results['p01'].append(p01)
            if not (np.isnan(effT_mK) or effT_mK < 0): 
                results['effT_mK'].append(effT_mK)   

    # 計算平均值和標準差
    p01_mean, p01_std = np.mean(results['p01']), np.std(results['p01'])
    effT_mean, effT_std = np.mean(results['effT_mK']), np.std(results['effT_mK'])
    
    # Convert results to a dictionary format suitable for plotting
    result_dict = {q: np.array(v) for q, v in results.items()}
    
    # 保存 histogram 圖片
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=save_folder)
    
    return folder_path, p01_mean, p01_std, effT_mean, effT_std

if __name__ == "__main__":
    QD_agent_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_25\DRKE#242_SumInfo.pkl"
    target_q = 'q0'
    base_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForRearrange\QDbackupIs1025"
    show_each_plot = False
    specific_folder_name="SS_30aveg" #"SS"
    # 初始化 Qmanager
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()

    # 用於存儲所有資料夾分析結果的列表
    analysis_results = []

    # 遍歷 base_folder 中的每個溫度資料夾
    for temp_folder in os.listdir(base_folder):
        temp_path = os.path.join(base_folder, temp_folder, specific_folder_name)
        
        if os.path.isdir(temp_path):
            # 呼叫分析函數並生成 histogram，並將結果添加到列表中
            analysis_results.append(analyze_and_plot(Qmanager, target_q, temp_path, show_each_plot, temp_path))
    
    # 統一列印所有資料夾的結果
    print("\nAll Folder Analysis Results:")
    for result in analysis_results:
        folder, p01_mean, p01_std, effT_mean, effT_std = result
        print(f"Folder: {folder}")
        print(f"p01 Mean: {p01_mean:.4f}, Std: {p01_std:.5f}")
        print(f"effT Mean: {effT_mean:.4f} mK, Std: {effT_std:.5f} mK\n")

    print("Analysis done! Check out the histogram in each SS folder.")

    

import numpy as np
import os
import re
from datetime import datetime

# Function to read dataset from a file
def read_dataset(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data = []
    for line in lines:
        date_str, time_str, value = line.strip().split(',')
        timestamp = datetime.strptime(f"{date_str} {time_str}", "%d-%m-%y %H:%M:%S")
        data.append((timestamp, float(value)))
    return data

# Function to filter data by timestamp range
def filter_data_by_time(data, start_time, end_time):
    return [value for timestamp, value in data if start_time <= timestamp <= end_time]

# Process the earliest and latest file timestamps
def extract_timestamp_from_filename(filename):
    match = re.search(r'H(\d{2})M(\d{2})S(\d{2})', filename)
    if match:
        hour, minute, second = map(int, match.groups())
        return hour, minute, second
    return None

# Main analysis function
def analyze_folder(folder_path, start_time, end_time):
    data_file = os.path.join(folder_path, "your_data_file.txt")  # Replace with actual file name
    data = read_dataset(data_file)

    # Filter data within the specified time range
    filtered_data = filter_data_by_time(data, start_time, end_time)
    
    # Convert to numpy array for analysis
    values = np.array(filtered_data) * 1e3  # Scale as per original script
    mean = np.mean(values)
    std_dev = np.std(values)

    print(f"Mean: {mean:.3f} +/- Std Dev: {std_dev:.4f}")
    return mean, std_dev

# Example of usage
folder_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForRearrange\QDbackupIs1025"
earliest_time = datetime(2024, 10, 27, 21, 32, 50)  # 27-10-24,21:32:50
latest_time = datetime(2024, 10, 27, 21, 43, 50)    # 27-10-24,21:43:50

# Analyze the folder data within the time range
analyze_folder(folder_path, earliest_time, latest_time)




