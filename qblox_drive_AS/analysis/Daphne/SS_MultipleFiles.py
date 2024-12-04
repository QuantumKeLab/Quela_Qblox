import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from qblox_drive_AS.support.UserFriend import *
from Modularize.support.QDmanager import QDmanager, Data_manager
from numpy import array, moveaxis, mean, std, median, arange
import numpy as np
from qblox_drive_AS.analysis.OneShotAna import a_OSdata_analPlot

def analyze_and_plot(Qmanager, target_q, folder_path, show_each_plot, save_folder):
    """針對特定資料夾中的 .nc 檔案進行分析並生成 histogram 圖"""
    results = {'p01': [], 'effT_mK': [], 'RO_fidelity_percentage':[], 'p10':[],'snr':[],'power_snr_dB':[],'dis':[],'sigma':[]}
    
    for file in os.listdir(folder_path):
        if file.endswith('.nc'):  # 確認檔案副檔名為 .nc
            nc_path = os.path.join(folder_path, file)
            p01, effT_mK, RO_fidelity_percentage, p10,snr,power_snr_dB,dis,sigma = a_OSdata_analPlot(Qmanager, target_q, nc_path, plot=show_each_plot)#
            
            # 檢查是否為有效值
            if not (np.isnan(p01) or p01 < 0):  
                results['p01'].append(p01)
            if not (np.isnan(effT_mK) or effT_mK < 0): 
                results['effT_mK'].append(effT_mK)   
            if not (np.isnan(RO_fidelity_percentage) or RO_fidelity_percentage < 0):  
                results['RO_fidelity_percentage'].append(RO_fidelity_percentage)
            if not (np.isnan(p10) or p10 < 0): 
                results['p10'].append(p10)   
            if not (np.isnan(snr) or snr < 0):  
                results['snr'].append(snr)
            if not (np.isnan(power_snr_dB) or power_snr_dB < 0): 
                results['power_snr_dB'].append(power_snr_dB)  
            if not (np.isnan(dis) or dis < 0): 
                results['dis'].append(dis)   
            if not (np.isnan(sigma) or sigma < 0):  
                results['sigma'].append(sigma)
            #It is "if", not "elif" (if you want to add more) 

    # 計算平均值和標準差
    p01_mean, p01_std = np.mean(results['p01']), np.std(results['p01'])
    effT_mean, effT_std = np.mean(results['effT_mK']), np.std(results['effT_mK'])
    RO_fidelity_percentage_mean, RO_fidelity_percentage_std = np.mean(results['RO_fidelity_percentage']), np.std(results['RO_fidelity_percentage'])
    p10_mean, p10_std = np.mean(results['p10']), np.std(results['p10'])
    snr_mean, snr_std = np.mean(results['snr']), np.std(results['snr'])
    power_snr_dB_mean, power_snr_dB_std = np.mean(results['power_snr_dB']), np.std(results['power_snr_dB'])
    dis_mean, dis_std = np.mean(results['dis']), np.std(results['dis'])
    sigma_mean, sigma_std = np.mean(results['sigma']), np.std(results['sigma'])

    # Convert results to a dictionary format suitable for plotting
    result_dict = {q: np.array(v) for q, v in results.items()}
    
    # 保存 histogram 圖片
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'p01', mode="pop", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'effT_mK', mode="ss", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'RO_fidelity_percentage', mode="rofdlty", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'p10', mode="p10", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'snr', mode="snr", save_fig=True, pic_folder=save_folder)
    # Data_manager().save_histo_pic(Qmanager, result_dict, 'power_snr_dB', mode="powersnr", save_fig=True, pic_folder=save_folder)
    
    return folder_path, p01_mean, p01_std, effT_mean, effT_std, RO_fidelity_percentage_mean, RO_fidelity_percentage_std, p10_mean, p10_std,snr_mean, snr_std, power_snr_dB_mean, power_snr_dB_std, sigma_mean, sigma_std, dis_mean, dis_std

if __name__ == "__main__":

    """Fill in"""
    QD_agent_path = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_26\DRKE#242_SumInfo.pkl"
    target_q = 'q0'
    main_folder = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\2024_10_26"
    show_each_plot = False
    specific_folder_name="SS_"#"SS" #"SS_30aveg"

    """Running"""
    # 初始化 Qmanager
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()

    # 用於存儲所有資料夾分析結果的列表
    analysis_results = []
    
    # 遍歷 main_folder 中的每個溫度資料夾
    for temp_folder in os.listdir(main_folder):
        temp_path = os.path.join(main_folder, temp_folder, specific_folder_name)
        
        if os.path.isdir(temp_path):
            # 呼叫分析函數並生成 histogram，並將結果添加到列表中
            analysis_results.append(analyze_and_plot(Qmanager, target_q, temp_path, show_each_plot, temp_path))
    # print('analysis_results:',analysis_results)
    # 統一列印所有資料夾的結果
    print("\nAll Folder Analysis Results:")
    for result in analysis_results:
        folder, p01_mean, p01_std, effT_mean, effT_std, RO_fidelity_percentage_mean, RO_fidelity_percentage_std, p10_mean, p10_std, snr_mean, snr_std, power_snr_dB_mean, power_snr_dB_std,sigma_mean, sigma_std, dis_mean, dis_std= result
        print(f"Folder: {folder}")
        print(f"p01 Mean: {p01_mean:.4f}, Std: {p01_std:.5f}")
        print(f"effT Mean: {effT_mean:.4f} mK, Std: {effT_std:.5f} mK")
        print(f"RO_fidelity Mean: {RO_fidelity_percentage_mean:.4f}% , Std: {RO_fidelity_percentage_std:.5f} %")
        print(f"p10 Mean: {p10_mean:.4f}%, Std: {p10_std:.5f}%")
        print(f"SNR Mean: {snr_mean:.4f}, Std: {snr_std:.5f}")
        print(f"Power SNR Mean: {power_snr_dB_mean:.4f}dB, Std: {power_snr_dB_std:.5f}dB")
        print(f"distance Mean: {dis_mean:.4f}, Std: {dis_std:.5f}")
        print(f"Sigma Mean: {sigma_mean:.4f}, Std: {sigma_std:.5f}\n")


    print("Analysis done! Check out the histograms in each SS folder.")


"""Clearer code but there is no unit shown in the terminal output result"""
# def analyze_and_plot(Qmanager, target_q, folder_path, show_each_plot, save_folder):
#     """針對特定資料夾中的 .nc 檔案進行分析並生成 histogram 圖"""
#     result_keys = ['p01', 'effT_mK', 'RO_fidelity_percentage', 'p10', 'snr', 'power_snr_dB']
#     results = {key: [] for key in result_keys}

#     # 分析每個 .nc 檔案
#     for file in os.listdir(folder_path):
#         if file.endswith('.nc'):
#             nc_path = os.path.join(folder_path, file)
#             analysis_results = a_OSdata_analPlot(Qmanager, target_q, nc_path, plot=show_each_plot)

#             # 遍歷結果並檢查有效值
#             for key, value in zip(result_keys, analysis_results):
#                 if not (np.isnan(value) or value < 0):
#                     results[key].append(value)

#     # 計算平均值和標準差
#     stats = {key: (np.mean(values), np.std(values)) for key, values in results.items()}

#     # 模式對應關係
#     mode_mapping = {
#         'p01': 'pop',
#         'effT_mK': 'ss',
#         'RO_fidelity_percentage': 'rofdlty',
#         'p10': 'p10',
#         'snr': 'snr',
#         'power_snr_dB': 'powersnr'
#     }

#     # 保存 histogram 圖片
#     data_manager = Data_manager()
#     for key in result_keys:
#         mode = mode_mapping[key]
#         data_manager.save_histo_pic(Qmanager, results, key, mode=mode, save_fig=True, pic_folder=save_folder)

#     return folder_path, stats



# if __name__ == "__main__":
#     """Fill in"""
#     QD_agent_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_26\DRKE#242_SumInfo.pkl"
#     target_q = 'q0'
#     main_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\TEST_FOLDER"
#     show_each_plot = False
#     specific_folder_name = "SS"  # "SS_30aveg"

#     """Running"""
#     # 初始化 Qmanager
#     Qmanager = QDmanager(QD_agent_path)
#     Qmanager.QD_loader()

#     # 用於存儲所有資料夾分析結果的列表
#     analysis_results = []

#     # 遍歷 main_folder 中的每個溫度資料夾
#     for temp_folder in os.listdir(main_folder):
#         temp_path = os.path.join(main_folder, temp_folder, specific_folder_name)

#         if os.path.isdir(temp_path):
#             # 呼叫分析函數並生成 histogram，並將結果添加到列表中
#             analysis_results.append(analyze_and_plot(Qmanager, target_q, temp_path, show_each_plot, temp_path))

#     # 統一列印所有資料夾的結果
#     print("\nAll Folder Analysis Results:")
#     for folder, stats in analysis_results:
#         print(f"Folder: {folder}")
#         for key, (mean, std) in stats.items():
#             print(f"{key} Mean: {mean:.4f}, Std: {std:.5f}")
#         print()

#     print("Analysis done! Check out the histograms in each SS folder.")