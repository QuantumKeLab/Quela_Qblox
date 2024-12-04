import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qcat.analysis.state_discrimination.readout_fidelity import GMMROFidelity, G1DROFidelity
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

def a_OSdata_analPlot(QD_agent:QDmanager, target_q:str, nc_path:str, plot:bool=True, pic_path:str='', save_pic:bool=False): # 
    folder = os.path.join(os.path.split(nc_path)[0],'OS_pic')
    if not os.path.exists(folder):
        os.mkdir(folder)
    if save_pic:
        pic_save_path = os.path.join(folder,os.path.split(nc_path)[1].split(".")[0]) if pic_path == '' else pic_path
    else:
        pic_save_path = None
    SS_ds = open_dataset(nc_path)
    ss_dict = Dataset.to_dict(SS_ds)
    # print(ss_dict)

    pe_I, pe_Q = ss_dict['data_vars']['e']['data']
    pg_I, pg_Q = ss_dict['data_vars']['g']['data']
   
    pgI_collection = [pg_I]
    pgQ_collection = [pg_Q]
    peI_collection = [pe_I]
    peQ_collection = [pe_Q]

    OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
    tarin_data, fit_arrays = OSdata_arranger(OS_data)
    # train GMM
    print(tarin_data)
    gmm2d_fidelity = GMMROFidelity()
    gmm2d_fidelity._import_data(tarin_data[0])
    gmm2d_fidelity._start_analysis()
    g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    print(gmm2d_fidelity.discriminator.cluster_model.means_)
    print(gmm2d_fidelity.label_map.label_assign.result)
    p00 = g1d_fidelity.g1d_dist[0][0][0]
    p01 = g1d_fidelity.g1d_dist[0][0][1]
    p10 = g1d_fidelity.g1d_dist[1][0][0]#
    p11 = g1d_fidelity.g1d_dist[1][0][1]
    p01=np.abs(p01)
    effT_mK = np.abs(p01_to_Teff(p01, transi_freq)*1000)#p01_to_Teff(p01, transi_freq)*1000
    RO_fidelity_percentage = (p00+p11)*100/2
    p10_precentage=p10*100
    snr = g1d_fidelity.discriminator.snr
    power_snr_dB=np.log10(snr)*20
    dis = g1d_fidelity.discriminator.signal
    sigma = np.mean(g1d_fidelity.discriminator.noise)
    # if plot:
    #     plot_readout_fidelity( tarin_data[0], gmm2d_fidelity, g1d_fidelity,transi_freq,pic_save_path, detail_output:bool=True)
    #     plt.close()
   
    if plot:
        # 新增變數賦值行，設定 detail_output=True
        fig, p01, effective_temp_mK, power_snr_dB,dis,sigma = plot_readout_fidelity(
            tarin_data[0], gmm2d_fidelity, g1d_fidelity, transi_freq, pic_save_path, plot=True, detail_output=True
        )
        
        # 現在你有了 `power_snr_dB` 和 `snr` 相關的值，可以進一步分析或保存
        print(f"Power SNR (dB): {power_snr_dB}")
        
        plt.close()

    print('|p01|',p01)
    print('effT_mK',effT_mK)
    print('p10',p10)
    print("snr",snr)
    print("power_snr_dB",power_snr_dB)
    print("distance",dis)
    print("sigma",sigma)
    return p01, effT_mK, RO_fidelity_percentage, p10_precentage,snr,power_snr_dB,dis,sigma

def a_OSdata_correlation_analPlot(nc_path:str, plot:bool=True, pic_path:str='', save_pic:bool=False): # 
    folder = os.path.join(os.path.split(nc_path)[0],'OS_pic')
    if not os.path.exists(folder):
        os.mkdir(folder)
    if save_pic:
        pic_save_path = os.path.join(folder,os.path.split(nc_path)[1].split(".")[0]) if pic_path == '' else pic_path
    else:
        pic_save_path = None
    SS_ds = open_dataset(nc_path)
    ss_dict = Dataset.to_dict(SS_ds)
    # print(ss_dict)
    
    # 讀取基態（g state）的 I 和 Q 數據
    pg_I = ss_dict['data_vars']['g']['data'][0]  # I 數據
    pg_Q = ss_dict['data_vars']['g']['data'][1]  # Q 數據

    # 將數據轉換為 numpy array 以便進行索引操作
    pg_I = np.array(pg_I)
    pg_Q = np.array(pg_Q)
        
    # 列出文件中的所有變數
    print(SS_ds.variables.keys())

    # 取得變數 'e' 和 'g'
    g_var = SS_ds.variables['g'][:]

    # 列出變數 'e' 和 'g' 的形狀以確認
    
    print(f"Shape of 'g': {g_var.shape}")
    # 確保數據總長度是偶數，以便進行分割
    assert len(pg_I) % 2 == 0, "數據總長度必須為偶數"

    # 使用 numpy 索引取奇數和偶數位置的數據
    I_readout_1 = pg_I[0::2]  # 從0開始每隔一個取一次
    I_readout_2 = pg_I[1::2] 
    
    Q_readout_1 = pg_Q[0::2]  
    Q_readout_2 = pg_Q[1::2] 

     # 將數據轉換為毫伏（選擇性）
    I_readout_1, I_readout_2 = array(I_readout_1) * 1000, array(I_readout_2) * 1000
    Q_readout_1, Q_readout_2 = array(Q_readout_1) * 1000, array(Q_readout_2) * 1000
    
    # I_readout_2_ave= np.mean(I_readout_2)
    # I_readout_1_ave= np.mean(I_readout_1)
    # Q_readout_2_ave= np.mean(Q_readout_2)
    # Q_readout_1_ave= np.mean(Q_readout_1)
    # 繪製 IQ 平面
    if plot:
        plt.figure(figsize=(6, 6))
        
        'Plot the scatter points for two readouts'
        plt.scatter(I_readout_1, Q_readout_1, c='blue', label="Readout 1", alpha=0.5, s=0.5)
        plt.scatter(I_readout_2, Q_readout_2, c='red', label="Readout 2", alpha=0.5,s=0.5)
        plt.axis('equal')

        # Label axes
        plt.xlabel('I (mV)')
        plt.ylabel('Q (mV)')
        plt.title('IQ Plane - Two Readouts')

        'Center=(0,0) and shows four quadrants'
        # # Set the axis limits symmetrically so that (0, 0) is at the center
        # max_lim = max(np.max(np.abs(I_readout_1)), np.max(np.abs(Q_readout_1)),
        #             np.max(np.abs(I_readout_2)), np.max(np.abs(Q_readout_2)))
        # plt.xlim(-max_lim, max_lim)
        # plt.ylim(-max_lim, max_lim)
        
        # # Add lines for the X and Y axes to show the quadrants
        # plt.axhline(0, color='black', linewidth=1.5)  # Horizontal line at y=0
        # plt.axvline(0, color='black', linewidth=1.5)  # Vertical line at x=0

        'Change ticks'
        # plt.tick_params(axis='both',which='major',labelsize=14)
        # x_major_locator=MultipleLocator(0.25)
        # #把x轴的刻度间隔设置为1,并存在变量里
        # y_major_locator=MultipleLocator(0.25)
        # #把y轴的刻度间隔设置为10,并存在变量里
        # ax=plt.gca()
        # #ax为两条坐标轴的实例
        # ax.xaxis.set_major_locator(x_major_locator)
        # #把x轴的主刻度设置为1的倍数
        # ax.yaxis.set_major_locator(y_major_locator)
        
        'Cycle'
        # 定義循環數
        # shots=10000
        # shots_per_cycle = 5
        # num_cycles = shots // shots_per_cycle  # 2000

        # 定義顏色列表，使用 5 種不同的顏色
        # colors = ['red','green', 'blue', 'orange' ,'purple']#,  

        # 創建一個圖形
        # plt.figure(figsize=(10, 6))

        ## 繪製數據
        # for cycle in range(num_cycles):
        #     for i in range(shots_per_cycle):
        #         idx = cycle * shots_per_cycle + i  # 計算出當前 shot 的索引
        #         plt.scatter(I_readout_1[idx], Q_readout_1[idx], color=colors[i], label=f'Shot {i+1}' if cycle == 0 else "")  # 只在第一個循環中顯示 label
        # for i in range(shots):
        #     plt.scatter(I_readout_1[i],Q_readout_1[i],color=colors[i],label=f'Shot{i+1}')
        
        # 添加標籤和標題
        # plt.xlim(106,108)
        # plt.ylim(-29,-25)
        # plt.xlabel('I Readout')
        # plt.ylabel('Q Readout')
        # plt.title('IQ Plane for Different Shots in a Cycle')

        # # 添加圖例（僅顯示一次每個顏色的標籤）
        # plt.legend(loc='best')

        'Plot'
        # Add grid and legend
        # plt.legend()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

        
        if save_pic:
            plt.savefig(f'{pic_save_path}_IQ_plane.png', bbox_inches='tight')
        
        plt.show()

    return I_readout_1, I_readout_2, Q_readout_1, Q_readout_2

def a_OSdata_analPlot_3G(QD_agent:QDmanager, target_q:str, nc_path:str, plot:bool=True, pic_path:str='', save_pic:bool=False): # 
    folder = os.path.join(os.path.split(nc_path)[0],'OS_pic')
    if not os.path.exists(folder):
        os.mkdir(folder)
    if save_pic:
        pic_save_path = os.path.join(folder,os.path.split(nc_path)[1].split(".")[0]) if pic_path == '' else pic_path
    else:
        pic_save_path = None
    SS_ds = open_dataset(nc_path)
    ss_dict = Dataset.to_dict(SS_ds)
    # print(ss_dict)
    pf_I, pf_Q = ss_dict['data_vars']['f']['data']
    pe_I, pe_Q = ss_dict['data_vars']['e']['data']
    pg_I, pg_Q = ss_dict['data_vars']['g']['data']

    pgI_collection = [pg_I]
    pgQ_collection = [pg_Q]
    peI_collection = [pe_I]
    peQ_collection = [pe_Q]
    pfI_collection = [pf_I]
    pfQ_collection = [pf_Q]

    OS_data = 1000*array([[pgI_collection,peI_collection, pfI_collection],[pgQ_collection,peQ_collection, pfQ_collection]]) # can train or predict 2*2*histo_counts*shot
    tarin_data, fit_arrays = OSdata_arranger(OS_data)
    # train GMM
    # print(tarin_data)
    gmm2d_fidelity = GMMROFidelity()
    gmm2d_fidelity._import_data(tarin_data[0])
    gmm2d_fidelity._start_analysis()
    g1d_fidelity = gmm2d_fidelity.export_G1DROFidelity()
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    print(gmm2d_fidelity.label_map.label_assign.result)
    p00 = g1d_fidelity.g1d_dist[0][0][0]
    p01 = g1d_fidelity.g1d_dist[0][0][1]
    p11 = g1d_fidelity.g1d_dist[1][0][1]
    effT_mK = np.abs(p01_to_Teff(p01, transi_freq)*1000)#p01_to_Teff(p01, transi_freq)*1000
    RO_fidelity_percentage = (p00+p11)*100/2
    if plot:
        plot_readout_fidelity( tarin_data[0], gmm2d_fidelity, g1d_fidelity,transi_freq,pic_save_path)
        plt.close()

    return p01, effT_mK, RO_fidelity_percentage

def share_model_OSana(QD_agent:QDmanager,target_q:str,folder_path:str,pic_save:bool=True):
    transi_freq = QD_agent.quantum_device.get_element(target_q).clock_freqs.f01()
    files = [name for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split("_")[1].split("(")[0]=="SingleShot")]
    files = [os.path.join(folder_path,name) for name in sort_files(files)][:21]
    pop_rec, efft_rec = [], []
    if pic_save:
        pic_folder = os.path.join(folder_path,"OS_detail_pic")
        if not os.path.exists(pic_folder):
            os.mkdir(pic_folder)
    else:
        pic_folder = None

    for idx, file in enumerate(files):
        SS_ds = open_dataset(file)
        ss_dict = Dataset.to_dict(SS_ds)
        # print(ss_dict)
        pe_I, pe_Q = ss_dict['data_vars']['e']['data']
        pg_I, pg_Q = ss_dict['data_vars']['g']['data']
        pgI_collection = [pg_I]
        pgQ_collection = [pg_Q]
        peI_collection = [pe_I]
        peQ_collection = [pe_Q]

        OS_data = 1000*array([[pgI_collection,peI_collection],[pgQ_collection,peQ_collection]]) # can train or predict 2*2*histo_counts*shot
        tarin_data, fit_arrays = OSdata_arranger(OS_data)
        # train GMM
        if idx == 0:
            gmm2d_fidelity = GMMROFidelity()
            gmm2d_fidelity._import_data(tarin_data[0])
            gmm2d_fidelity._start_analysis()

        gmm2d_fidelity.discriminator._import_data( fit_arrays[0] )
        gmm2d_fidelity.discriminator._start_analysis()


    return pop_rec, efft_rec





if __name__ == "__main__":

    
    QD_agent_path=r"C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_12_4\DRKE#242_SumInfo.pkl"
    Qmanager = QDmanager(QD_agent_path)
    Qmanager.QD_loader()
    target_q='q0'

    "For single file"
    nc_path=r"C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_12_4\DRKEq0_SingleShot(0)_H21M53S12.nc"
    a_OSdata_analPlot(Qmanager,target_q, nc_path)
    # a_OSdata_correlation_analPlot(nc_path)
    
