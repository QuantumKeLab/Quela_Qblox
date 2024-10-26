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

def generate_histogram(data, title, save_path=None):
    plt.figure()
    plt.hist(data, bins=30, alpha=0.7, color='b')
    plt.title(title)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    if save_path:
        plt.savefig(save_path)
    plt.close()

def share_model_OSana(QD_agent, target_q, folder_path, pic_save=True):
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    files = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path, name)) and name.split("_")[1].split("(")[0] == "SingleShot")]
    files = [os.path.join(folder_path, name) for name in sort_files(files)][:21]
    pop_rec, efft_rec = [], []
    
    if pic_save:
        pic_folder = os.path.join(folder_path, "OS_detail_pic")
        if not os.path.exists(pic_folder):
            os.mkdir(pic_folder)
    else:
        pic_folder = None

    # 初始化 GMM 2D 模型
    gmm2d_fidelity = GMMROFidelity()
    
    for idx, file in enumerate(files):
        SS_ds = open_dataset(file)
        ss_dict = SS_ds.to_dict()
        
        pe_I, pe_Q = ss_dict['data_vars']['e']['data']
        pg_I, pg_Q = ss_dict['data_vars']['g']['data']
        
        OS_data = 1000 * array([[pg_I, pe_I], [pg_Q, pe_Q]])
        train_data, fit_arrays = OSdata_arranger(OS_data)
        
        # 進行 GMM 訓練與分析
        if idx == 0:
            gmm2d_fidelity._import_data(train_data[0])
            gmm2d_fidelity._start_analysis()
        
        gmm2d_fidelity.discriminator._import_data(fit_arrays[0])
        gmm2d_fidelity.discriminator._start_analysis()
        
        # 將結果保存
        pop_rec.append(gmm2d_fidelity.pop_value)
        efft_rec.append(gmm2d_fidelity.efficiency_value)
    
    # 生成並保存 histogram 圖
    if pic_save and pic_folder:
        pop_hist_path = os.path.join(pic_folder, f"{target_q}_population_histogram.png")
        efft_hist_path = os.path.join(pic_folder, f"{target_q}_efficiency_histogram.png")
        
        generate_histogram(pop_rec, f"{target_q} Population Histogram", pop_hist_path)
        generate_histogram(efft_rec, f"{target_q} Efficiency Histogram", efft_hist_path)
    
    return pop_rec, efft_rec

if __name__ == "__main__":
    QD_agent_path = r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_10_25\DRKE#242_SumInfo.pkl"
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()
    
    target_q = 'q0'
    folder_path = r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_25\15mK_IntegrationTime_1250ns\SS"
    
    efft_rec = []
    thermal_pop = []
    
    info = share_model_OSana(Qmanager, target_q, folder_path)
    efft_rec.append(info[1])
    thermal_pop.append(info[0])
    
    highlight_print(f"{target_q}: {round(median(array(efft_rec)), 2)} +/- {round(std(array(efft_rec)), 3)} mK")
    Data_manager().save_histo_pic(Qmanager, efft_rec, target_q, mode="ss")
    Data_manager().save_histo_pic(Qmanager, thermal_pop, target_q, mode="pop")
