import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from xarray import Dataset, open_dataset

from qblox_drive_AS.support.UserFriend import *
from qblox_drive_AS.support.QDmanager import QDmanager, Data_manager
from numpy import array, moveaxis, mean, std, median, arange
import matplotlib.pyplot as plt

import numpy as np
from matplotlib.ticker import MultipleLocator
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set
from qblox_drive_AS.analysis.OneShotAna import a_OSdata_analPlot
import re
from datetime import datetime, timedelta


def extract_datetime_from_filename(filename):
    """從檔名中提取時間，假設格式為 H18M25S41 並轉換為 datetime 對象"""
    time_str = re.search(r'H(\d{2})M(\d{2})S(\d{2})', filename)
    if time_str:
        hours, minutes, seconds = map(int, time_str.groups())
        today = datetime.now().replace(hour=hours, minute=minutes, second=seconds, microsecond=0)
        return today
    return None

def analyze_and_plot(Qmanager, target_q, folder_path, show_each_plot, save_folder, time_mode="relative"):
    """針對特定資料夾中的 .nc 檔案進行分析並生成 histogram 和隨時間變化的圖"""
    results = {'p01': [], 'effT_mK': [], 'time': [], 'real_time': []}
    
    for file in os.listdir(folder_path):
        if file.endswith('.nc'):  # 確認檔案副檔名為 .nc
            nc_path = os.path.join(folder_path, file)
            p01, effT_mK, RO_fidelity_percentage, p10,snr,power_snr_dB= a_OSdata_analPlot(Qmanager, target_q, nc_path, plot=show_each_plot)

            # 提取檔案中的時間資訊
            file_datetime = extract_datetime_from_filename(file)
            if file_datetime:
                results['real_time'].append(file_datetime)

            # 檢查 p01 和 effT_mK 是否為有效值
            if not (np.isnan(p01) or p01 < 0):  
                results['p01'].append(p01)
            if not (np.isnan(effT_mK) or effT_mK < 0): 
                results['effT_mK'].append(effT_mK)

    # 根據最早的時間點來計算相對時間
    if results['real_time']:
        start_time = min(results['real_time'])
        results['time'] = [(t - start_time).total_seconds() for t in results['real_time']]

    # 確保所有數據長度一致
    min_length = min(len(results['p01']), len(results['effT_mK']), len(results['time']))
    results = {key: values[:min_length] for key, values in results.items()}

    # 計算平均值和標準差
    p01_mean, p01_std = np.mean(results['p01']), np.std(results['p01'])
    effT_mean, effT_std = np.mean(results['effT_mK']), np.std(results['effT_mK'])

    # 保存 histogram 圖片
    result_dict = {q: np.array(v) for q, v in results.items()}
    Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=save_folder)
    Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=save_folder)

    # 畫出 effT 和 p01 隨時間變化的圖
    if results['time']:
        x_axis = results['time']
        x_label = "Time (seconds from start)"

        # 建立圖表，顯示所有數據
        plt.figure(figsize=(12, 6))

        # p01 vs time
        plt.plot(x_axis, results['p01'], '-o', label='p01')

        # effT vs time
        plt.plot(x_axis, results['effT_mK'], '-o', label='effT (mK)')

        plt.xlabel(x_label)
        plt.ylabel("Values")
        plt.title("p01 and effT over Time")
        plt.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, f"effT_p01_time_variation_{time_mode}.png"))
        plt.show()

    return folder_path, p01_mean, p01_std, effT_mean, effT_std

if __name__ == "__main__":
    QD_agent_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_24\DRKE#242_SumInfo.pkl"
    target_q = 'q0'
    base_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupCouldBe1024\SS_overnight"
    show_each_plot = False
    time_mode = "relative"  # "real" 或 "relative"

    # 初始化 Qmanager
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()

    # 用於存儲所有資料夾分析結果的列表
    analysis_results = []

    # 遍歷 base_folder 中的每個溫度資料夾
    for temp_folder in os.listdir(base_folder):
        temp_path = os.path.join(base_folder, temp_folder, "SS")
        
        if os.path.isdir(temp_path):
            # 呼叫分析函數並生成 histogram，並將結果添加到列表中
            analysis_results.append(analyze_and_plot(Qmanager, target_q, temp_path, show_each_plot, temp_path, time_mode))
    
    # 統一列印所有資料夾的結果
    print("\nAll Folder Analysis Results:")
    for result in analysis_results:
        folder, p01_mean, p01_std, effT_mean, effT_std = result
        print(f"Folder: {folder}")
        print(f"p01 Mean: {p01_mean:.4f}, Std: {p01_std:.5f}")
        print(f"effT Mean: {effT_mean:.4f} mK, Std: {effT_std:.5f} mK\n")

    print("Analysis done! Check out the histogram and time variation plots in each SS folder.")

