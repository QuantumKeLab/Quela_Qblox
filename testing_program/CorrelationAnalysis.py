from xarray import open_dataset
from xarray import Dataset
from numpy import array
import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd

"""Get ~Vg and ~Ve from RunI and RunII, respectively"""
def ave_Vg_and_ave_Ve(ave_V_nc_path:str):
    ave_V_nc_path = r'C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_1\DRKEq0_SingleShot(0)_H22M15S28.nc'
    dataset = nc.Dataset(ave_V_nc_path, 'r')

    # 列出文件中的所有變數
    print(dataset.variables.keys())

    # 或者詳細列出每個變數的信息
    for var_name in dataset.variables:
        print(f"Variable: {var_name}, Shape: {dataset.variables[var_name].shape}, Dimensions: {dataset.variables[var_name].dimensions}")


    # 取得變數 'e' 和 'g'
    e_var = dataset.variables['e'][:]
    g_var = dataset.variables['g'][:]

    # 列出變數 'e' 和 'g' 的形狀以確認
    print(f"Shape of 'e': {e_var.shape}")
    print(f"Shape of 'g': {g_var.shape}")

    # # 列印變數 'e' 和 'g' 的原始資料（平均前的數據）
    # print(f"Data of 'e' before averaging:\n{e_var[:,:10]}")
    # print(f"Data of 'g' before averaging:\n{g_var[:,:10]}")
    e_I0= e_var[0,:]
    g_I0= g_var[0,:]
    # 對第一個維度 ('I') 進行平均， axis=0 表示沿著第 0 個維度進行操作
    Ve_mean_I = np.mean(e_I0) #|e>的I
    Vg_mean_I = np.mean(g_I0) #|g>的I

    # 列印出平均後的結果
    print(f"Mean of 'e' along dimension 'I': {Ve_mean_I}")
    print(f"Mean of 'g' along dimension 'I': {Vg_mean_I}")

    # 關閉文件
    dataset.close()
    return Ve_mean_I, Vg_mean_I

"""Get Vk_0 and Vk_tau from RunI """
def Vk_0_and_Vk_tau(Vk_nc_file:str):
    # Open the NetCDF file
    Vk_nc_file = r'C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_1\DRKEq0_SingleShot(0)_H22M16S18.nc'
    ds = open_dataset(Vk_nc_file)
    ss_dict = Dataset.to_dict(ds)
    print(ds.data_vars)

    print(array(ss_dict['data_vars']['g']['data']).shape)

    # 讀取基態（g state）的 I 和 Q 數據
    pg_I = ss_dict['data_vars']['g']['data'][0]  # I 數據
    pg_Q = ss_dict['data_vars']['g']['data'][1]  # Q 數據

    # 將數據轉換為 numpy array 以便進行索引操作
    pg_I = np.array(pg_I)
    pg_Q = np.array(pg_Q)

    # 確保數據總長度是偶數，以便進行分割
    assert len(pg_I) % 2 == 0, "數據總長度必須為偶數"

    # 使用 numpy 索引取奇數和偶數位置的數據
    I_readout_2 = pg_I[0::2]  # 取奇數位 (第一次 readout)
    I_readout_1 = pg_I[1::2] # 取偶數位 (第二次 readout)

    Q_readout_2 = pg_Q[0::2]  # 取奇數位 (第一次 readout)
    Q_readout_1 = pg_Q[1::2] # 取偶數位 (第二次 readout)

    # 將它們轉換成 list
    I_readout_2_list =I_readout_2.tolist()
    I_readout_1_list =I_readout_1.tolist()
    Q_readout_2_list =Q_readout_2.tolist()
    Q_readout_1_list =Q_readout_1.tolist()


    # 列印結果 
    # print("First readout of I: ", I_readout_1_list)
    # print("******")
    # print("Second readout of I: ", I_readout_2_list)
    # print("******")
    return I_readout_2_list, I_readout_1_list, Q_readout_2_list, Q_readout_1_list

"""Get normalized_g_1, Pe and Teff"""
def g_1(I_readout_2_list, I_readout_1_list,Ve_mean_I, Vg_mean_I):
    # 假設 Vk_0 和 Vk_tau 是你從前面的數據中提取的兩個 list
    I_readout_2_list, I_readout_1_list=Vk_0_and_Vk_tau(Vk_nc_file='')
    Vk_0 = I_readout_1_list  # 替換為你的實際 list
    Vk_tau = I_readout_2_list  # 替換為你的實際 list

    # 假設你已經有 V̄e 和 V̄g，這些是常數
    Ve_mean_I, Vg_mean_I =ave_Vg_and_ave_Ve(ave_V_nc_path='')
    V_e_avg = Ve_mean_I  
    V_g_avg = Vg_mean_I  

    # 計算歸一化的 V̄e 和 V̄g 差異
    V_diff = V_e_avg - V_g_avg

    # 確保 Vk_0 和 Vk_tau 的長度相同
    N = len(Vk_0) # should be equal to shot number

    # 用於存放累積的結果
    sum_result = 0

    # 逐步計算公式
    for k in range(N):
        term_0 = (Vk_0[k] - V_g_avg) / V_diff  # Re[(Vk(0) - V̄g) / (V̄e - V̄g)]
        term_tau = (Vk_tau[k] - V_g_avg) / V_diff  # Re[(Vk(τ) - V̄g) / (V̄e - V̄g)]
        
        # 對每一項取實部，並計算積的累積
        sum_result += np.real(term_0) * np.real(term_tau)

    # 最終結果除以 N
    normalized_g_1 = sum_result / N

    # 列印結果
    print("Normalized g^(1)= ", normalized_g_1)
    return normalized_g_1, Vk_0, Vk_tau, V_e_avg, V_g_avg

# """"""
if __name__ == '__main__':
    normalized_g_1= g_1(normalized_g_1)
    # N=shot_num
    f01=4e9
    # normalized_g_1=g_1/N  #??
    P_e= 1/2-1/(2*np.sqrt(1+4*normalized_g_1))
    qubit_frequency=2*np.pi*f01
    h_bar = 1.054571800*1e-34
    k_B = 1.38e-23
    T_eff=-h_bar*qubit_frequency/(k_B*np.log(P_e/(1-P_e)))*1000

    print(f"effective temperature = {T_eff}mK")