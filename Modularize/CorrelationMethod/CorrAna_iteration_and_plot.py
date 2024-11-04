import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import re
from datetime import datetime
from lmfit import Model

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

def T1_func(D, A, T1, offset):
    return A * np.exp(-D / T1) + offset

def T1_fit_analysis(data: np.ndarray, freeDu: np.ndarray, T1_guess: float = 10 * 1e-6, return_error: bool = False):
    offset_guess = data[-1]
    a_guess = np.max(data) - offset_guess if data[-1] > offset_guess else np.min(data) - offset_guess
    T1_model = Model(T1_func)
    result = T1_model.fit(data, D=freeDu, A=a_guess, T1=T1_guess, offset=offset_guess)

    A_fit = result.best_values['A']
    T1_fit = result.best_values['T1']
    offset_fit = result.best_values['offset']
    para_fit = np.linspace(freeDu.min(), freeDu.max(), 50 * len(data))
    fitting = T1_func(para_fit, A_fit, T1_fit, offset_fit)
    
    if not return_error:
        return xr.Dataset(
            data_vars=dict(
                data=(['freeDu'], data),
                fitting=(['para_fit'], fitting)
            ),
            coords=dict(
                freeDu=(['freeDu'], freeDu),
                para_fit=(['para_fit'], para_fit)
            ),
            attrs=dict(
                exper="T1",
                T1_fit=T1_fit
            )
        )
    else:
        fit_error = float(result.covar[1][1]) * 1e6
        return xr.Dataset(
            data_vars=dict(
                data=(['freeDu'], data),
                fitting=(['para_fit'], fitting)
            ),
            coords=dict(
                freeDu=(['freeDu'], freeDu),
                para_fit=(['para_fit'], para_fit)
            ),
            attrs=dict(
                exper="T1",
                T1_fit=T1_fit
            )
        ), fit_error

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
    
    # Step 3: Fit T1 and plot results
    T1_dataset = T1_fit_analysis(g1_values, correlate_delays)

    plt.figure(figsize=(10, 6))
    plt.plot(correlate_delays, g1_values, 'o', label="g^(1) Data")
    plt.plot(T1_dataset['para_fit'], T1_dataset['fitting'], '-', label="T1 Fit")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("g^(1)")
    plt.legend()
    plt.grid(True)
    
    for delay, teff in zip(correlate_delays, teff_values[1:]):
        print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")

    plt.show()
    return


if __name__ == '__main__':
    ave_Ve_nc_file = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026\40mK\SS\DRKEq0_SingleShot(16)_H14M20S31.nc"
    Vk_directory = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026\40mK\SS_corre"
    correlate_delays = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40]
    f01 = 4.31791e9
    correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays)

