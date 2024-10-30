import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

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

def correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays):
    Ve_mean_I = ave_Ve(ave_Ve_nc_file)
    
    # Step 1: Group and sort .nc files
    Vk_files_grouped = {0: [], 1: [], 2: []}
    for file in os.listdir(Vk_directory):
        if file.endswith('.nc'):
            index = int(file.split("_SingleShot(")[1].split(")_")[0])
            Vk_files_grouped[index].append(os.path.join(Vk_directory, file))

    for key in Vk_files_grouped:
        Vk_files_grouped[key].sort()  # Sort files in each group by time

    # Step 2: Calculate average g1 for each group
    g1_values_per_group = []
    teff_values_per_group = []
    
    for index, files in Vk_files_grouped.items():
        g1_values = []
        teff_values = []
        ave_Vg_nc_file = files[0]  # First file in each group
        Vg_mean_I = ave_Vg(ave_Vg_nc_file)

        for Vk_nc_file in files:
            g1, teff = calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I, f01)
            g1_values.append(g1)
            teff_values.append(teff)

        # Average g1 for each group
        g1_values_per_group.append(np.mean(g1_values))
        teff_values_per_group.append(np.mean(teff_values))

    # Step 3: Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(correlate_delays, g1_values_per_group, marker='o', color='b', label="Average g^(1) per group")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("Average g^(1)")
    plt.legend()
    plt.grid(True)
    
    for delay, teff in zip(correlate_delays, teff_values_per_group):
        print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")

    plt.show()
    return

if __name__ == '__main__':
    ave_Ve_nc_file = r"C:\Users\Ke Lab\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026\55mK\SS\DRKEq0_SingleShot(31)_H21M23S46.nc"
    Vk_directory = r"C:\Users\Ke Lab\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026\55mK\SS_corre"
    correlate_delays = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40]
    f01 = 4.31791e9
    correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays)
