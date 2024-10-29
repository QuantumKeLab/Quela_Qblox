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


######可用但是跑長時間的，圖只出現前面幾筆####
from datetime import datetime, timedelta



def extract_time_from_filename(filename):
    """從檔名中提取時間並轉換為時間字串和秒數，格式假設為 H18M25S41"""
    time_str = re.search(r'H(\d{2})M(\d{2})S(\d{2})', filename)
    if time_str:
        hours, minutes, seconds = map(int, time_str.groups())
        file_time = timedelta(hours=hours, minutes=minutes, seconds=seconds)
        total_seconds = file_time.total_seconds()
        return file_time, total_seconds
    return None, None

def analyze_and_plot(Qmanager, target_q, folder_path, show_each_plot, save_folder, time_mode="relative"):
    """針對特定資料夾中的 .nc 檔案進行分析並生成 histogram 和隨時間變化的圖"""
    results = {'p01': [], 'effT_mK': [], 'time': [], 'real_time': []}
    
    for file in os.listdir(folder_path):
        if file.endswith('.nc'):  # 確認檔案副檔名為 .nc
            nc_path = os.path.join(folder_path, file)
            p01, effT_mK, _ = a_OSdata_analPlot(Qmanager, target_q, nc_path, plot=show_each_plot)

            # 提取檔案中的時間資訊
            file_time, time_in_seconds = extract_time_from_filename(file)
            if time_in_seconds is not None:
                results['time'].append(time_in_seconds)
                results['real_time'].append(file_time)

            # 檢查 p01 和 effT_mK 是否為有效值
            if not (np.isnan(p01) or p01 < 0):  
                results['p01'].append(p01)
            if not (np.isnan(effT_mK) or effT_mK < 0): 
                results['effT_mK'].append(effT_mK)   

    # 計算平均值和標準差
    p01_mean, p01_std = np.mean(results['p01']), np.std(results['p01'])
    effT_mean, effT_std = np.mean(results['effT_mK']), np.std(results['effT_mK'])
    
    # Convert results to a dictionary format suitable for plotting histograms
    result_dict = {q: np.array(v) for q, v in results.items()}
    
    # 保存 histogram 圖片
    Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=save_folder)
    Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=save_folder)

    # 畫出 effT 和 p01 隨時間變化的圖
    if results['time']:
        times_sorted, p01_sorted, effT_sorted = zip(*sorted(zip(results['time'], results['p01'], results['effT_mK'])))
        real_times_sorted = [time.total_seconds() for time in sorted(results['real_time'])]

        # 設定 x 軸顯示模式
        if time_mode == "real":
            x_axis = real_times_sorted
            x_ticks = [str(timedelta(seconds=int(t))) for t in x_axis]  # 格式化為 HH:MM:SS
            x_label = "Time (HH:MM:SS)"
        else:
            x_axis = np.array(times_sorted) - times_sorted[0]
            x_ticks = x_axis
            x_label = "Time (seconds from start)"

        plt.figure(figsize=(10, 5))
        
        # p01 vs time
        plt.subplot(1, 2, 1)
        plt.plot(x_axis, p01_sorted, '-o')
        plt.xlabel(x_label)
        plt.ylabel("p01")
        plt.title("p01 over Time")
        if time_mode == "real":
            plt.xticks(ticks=x_axis, labels=x_ticks, rotation=45)
        
        # effT vs time
        plt.subplot(1, 2, 2)
        plt.plot(x_axis, effT_sorted, '-o')
        plt.xlabel(x_label)
        plt.ylabel("effT (mK)")
        plt.title("effT over Time")
        if time_mode == "real":
            plt.xticks(ticks=x_axis[::max(len(x_axis)//10, 1)], labels=x_ticks[::max(len(x_ticks)//10, 1)], rotation=45)

        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, f"effT_p01_time_variation_{time_mode}.png"))
        plt.show()
    
    return folder_path, p01_mean, p01_std, effT_mean, effT_std

if __name__ == "__main__":
    QD_agent_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_26\DRKE#242_SumInfo.pkl"
    target_q = 'q0'
    base_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026OrCouldBe1028"
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

