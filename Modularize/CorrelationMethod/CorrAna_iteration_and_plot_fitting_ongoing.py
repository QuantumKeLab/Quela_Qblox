import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import re
from datetime import datetime
from scipy.optimize import curve_fit

def ave_Ve(ave_V_nc_file: str):
    dataset = xr.open_dataset(ave_V_nc_file)
    e_var = dataset['e'].values
    e_I0 = e_var[0, :]
    Ve_mean_I = np.mean(e_I0)
    dataset.close()
    return Ve_mean_I

def ave_Vg(ave_Vg_nc_file: str):
    ds = xr.open_dataset(ave_Vg_nc_file)
    pg_I = ds['g'].values[0]
    I_readout_1 = pg_I[0::2]
    Vg_mean_I = np.mean(I_readout_1)
    ds.close()
    return Vg_mean_I

def Vk_0_and_Vk_tau(Vk_nc_file: str):
    ds = xr.open_dataset(Vk_nc_file)
    pg_I = ds['g'].values[0]
    I_readout_1 = pg_I[0::2]
    I_readout_2 = pg_I[1::2]
    ds.close()
    return I_readout_1.tolist(), I_readout_2.tolist()

def g_1(I_readout_1_list, I_readout_2_list, Ve_mean_I, Vg_mean_I):
    Vk_0 = I_readout_1_list
    Vk_tau = I_readout_2_list
    V_diff = Ve_mean_I - Vg_mean_I
    N = len(Vk_0)
    sum_result = sum((Vk_0[k] - Vg_mean_I) / V_diff * (Vk_tau[k] - Vg_mean_I) / V_diff for k in range(N))
    normalized_g_1 = sum_result / N
    return normalized_g_1

def calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I, f01):
    I_readout_1_list, I_readout_2_list = Vk_0_and_Vk_tau(Vk_nc_file)
    normalized_g_1 = g_1(I_readout_1_list, I_readout_2_list, Ve_mean_I, Vg_mean_I)
    P_e = (1 / 2) - 1 / (2 * np.sqrt(1 + 4 * normalized_g_1))
    qubit_frequency = 2 * np.pi * f01
    h_bar = 1.054571800e-34
    k_B = 1.381e-23
    T_eff = -h_bar * qubit_frequency / (k_B * np.log(P_e / (1 - P_e))) * 1000
    return normalized_g_1, T_eff

def t1(D, A, T1, offset):
    return A * np.exp(-D / T1) + offset

def correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays):
    Ve_mean_I = ave_Ve(ave_Ve_nc_file)
    
    # Step 1: Sort .nc files by time, only include index 0
    Vk_files = [os.path.join(Vk_directory, file) for file in os.listdir(Vk_directory) if file.endswith('.nc') and "_SingleShot(0)_" in file]
    Vk_files.sort(key=lambda x: datetime.strptime(re.search(r'_H(\d+)M(\d+)S(\d+)', x).group(0), '_H%HM%MS%S'))

    # Step 2: Calculate g1 and Teff for each file
    g1_values = []
    teff_values = []
    ave_Vg_nc_file = Vk_files[0]  # Use the first file for Vg_mean_I
    print(f"Using Vg file: {ave_Vg_nc_file}")
    Vg_mean_I = ave_Vg(ave_Vg_nc_file)

    for Vk_nc_file in Vk_files:
        print(f"Processing Vk file: {Vk_nc_file}")
        g1, teff = calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I, f01)
        g1_values.append(g1)
        teff_values.append(teff)

    # Exclude the first data point
    g1_values = np.array(g1_values[1:])
    correlate_delays = np.array(correlate_delays[1:len(g1_values) + 1])
    
    # Step 3: Fit T1 using curve_fit from scipy
# 調整 A 和 offset 的初始猜測值和邊界
    A_guess = g1_values[0] - g1_values[-1]
    T1_guess = np.mean(correlate_delays)
    offset_guess = g1_values[-1]

    # 增加 A 和 offset 的合理上下限範圍
    A_guess_upper = max(2 * (g1_values[0] - g1_values[-1]), 0.5 * (g1_values[0] - g1_values[-1]))
    A_guess_lower = min(-2 * abs(g1_values[0] - g1_values[-1]), 0.5 * (g1_values[0] - g1_values[-1]))  # 允許 A 為負
    C_guess_upper = max(0.1 * g1_values[-1], 2 * g1_values[-1])
    C_guess_lower = 0  # 假設 offset 最小不能小於零

    bounds = ((A_guess_lower, 0.1 * T1_guess, C_guess_lower), (A_guess_upper, 3 * T1_guess, C_guess_upper))
    p0 = (A_guess, T1_guess, offset_guess)
    ans, ans_error = curve_fit(t1, correlate_delays, g1_values, p0=p0, bounds=bounds)

    # Calculate effective temperature when g^(1) = 0
    T1_fit, offset_fit = ans[1], ans[2]
    g1_at_y0_delay = -T1_fit * np.log(-offset_fit / ans[0]) if ans[0] < 0 else None

    if g1_at_y0_delay is not None:
        print(f"g^(1) reaches 0 at delay: {g1_at_y0_delay:.2f} µs")
    else:
        print("g^(1) does not reach 0 within the fitted range.")

    # Step 4: Plot results
    plt.figure(figsize=(10, 6))
    plt.scatter(correlate_delays, g1_values, c='blue', label="g^(1) Data")
    plt.plot(correlate_delays, t1(correlate_delays, *ans), c='red', label="T1 Fit")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("g^(1)")
    plt.legend()
    plt.grid(True)
    
    for delay, teff in zip(correlate_delays, teff_values[1:]):
        print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")

    plt.show()
    return

if __name__ == '__main__':
    ave_Ve_nc_file = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026OrCouldBe1028\100mK\SS\DRKEq0_SingleShot(29)_H17M37S18.nc"
    Vk_directory = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026OrCouldBe1028\100mK\SS_corre"
    correlate_delays = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40]
    f01 = 4.31791e9
    correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays)
