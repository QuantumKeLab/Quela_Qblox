from xarray import open_dataset
from xarray import Dataset
from numpy import array
# import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd

"""Get ~Vg and ~Ve from RunI and RunII, respectively"""
def ave_Vg_and_ave_Ve(ave_V_nc_path:str):
    # dataset = nc.Dataset(ave_V_nc_path, 'r')
    dataset = open_dataset(ave_V_nc_path)
    
    # 列出文件中的所有變數
    print(dataset.variables.keys())

    # 取得變數 'e' 和 'g'
    e_var = dataset.variables['e'][:]
    g_var = dataset.variables['g'][:]

    # 列出變數 'e' 和 'g' 的形狀以確認
    print(f"Shape of 'e': {e_var.shape}")
    print(f"Shape of 'g': {g_var.shape}")

    # 取得基態和激發態的 I 方向數據並取平均
    e_I0 = e_var[0,:]
    g_I0 = g_var[0,:]
    Ve_mean_I = np.mean(e_I0)  # |e> 的 I 平均
    Vg_mean_I = np.mean(g_I0)  # |g> 的 I 平均

    print(f"Mean of 'e' along dimension 'I': {Ve_mean_I}")
    print(f"Mean of 'g' along dimension 'I': {Vg_mean_I}")

    # 關閉文件
    dataset.close()
    return Ve_mean_I, Vg_mean_I

"""Get Vk_0 and Vk_tau from RunI """
def Vk_0_and_Vk_tau(Vk_nc_file:str):
    ds = open_dataset(Vk_nc_file)
    ss_dict = Dataset.to_dict(ds)

    # 讀取基態（g state）的 I 數據
    pg_I = np.array(ss_dict['data_vars']['g']['data'][0])

    # 確保數據總長度是偶數，以便進行分割
    assert len(pg_I) % 2 == 0, "數據總長度必須為偶數"

    # 使用 numpy 索引取奇數和偶數位置的數據
    I_readout_1 = pg_I[0::2]  # 取奇數位 (第一次 readout)
    I_readout_2 = pg_I[1::2]  # 取偶數位 (第二次 readout)

    # 返回 list 格式
    return I_readout_1.tolist(), I_readout_2.tolist()

"""Get normalized_g_1, Pe and Teff"""
def g_1(I_readout_1_list, I_readout_2_list, Ve_mean_I, Vg_mean_I):
    Vk_0 = I_readout_1_list  # 使用前面的讀數
    Vk_tau = I_readout_2_list

    # 計算歸一化的 V̄e 和 V̄g 差異
    V_diff = Ve_mean_I - Vg_mean_I
    N = len(Vk_0)  # 數據點數量

    sum_result = 0
    for k in range(N):
        term_0 = (Vk_0[k] - Vg_mean_I) / V_diff
        term_tau = (Vk_tau[k] - Vg_mean_I) / V_diff
        sum_result += np.real(term_0) * np.real(term_tau)

    normalized_g_1 = sum_result / N
    print("Normalized g^(1)= ", normalized_g_1)
    return normalized_g_1

if __name__ == '__main__':
    # Step 1: Get Ve_mean_I and Vg_mean_I
    ave_V_nc_path = r'C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_9\cor\DR2q0_SingleShot(0)_H10M18S59.nc'
    Ve_mean_I, Vg_mean_I = ave_Vg_and_ave_Ve(ave_V_nc_path)
    # Ve_mean_I=5e-5
    # Vg_mean_I=5.5e-5

    # Step 2: Get I_readout_2_list and I_readout_1_list
    Vk_nc_file = r'C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_9\cor\DR2q0_SingleShot(0)_H10M26S42.nc'
    I_readout_1_list, I_readout_2_list = Vk_0_and_Vk_tau(Vk_nc_file)

    # Step 3: Calculate normalized g^(1)
    normalized_g_1 = g_1(I_readout_1_list, I_readout_2_list, Ve_mean_I, Vg_mean_I)

    # Step 4: Calculate P_e and T_eff
    f01 = 4e9  # qubit frequency
    P_e = 1 / 2 - 1 / (2 * np.sqrt(1 + 4 * normalized_g_1))

    qubit_frequency = 2 * np.pi * f01
    h_bar = 1.054571800e-34
    k_B = 1.38e-23

    T_eff = -h_bar * qubit_frequency / (k_B * np.log(P_e / (1 - P_e))) * 1000
    print(f"Effective temperature = {T_eff} mK")
