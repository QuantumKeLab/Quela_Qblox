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

def correlation_method(ave_Ve_nc_file, Vk_directory, f01, correlate_delays, plot_title):
    print(f"Using ave_Ve file: {ave_Ve_nc_file}")
    Ve_mean_I = ave_Ve(ave_Ve_nc_file)
    
    # Step 1: Sort .nc files by time, regardless of index
    Vk_files = [os.path.join(Vk_directory, file) for file in os.listdir(Vk_directory) if file.endswith('.nc')]
    Vk_files.sort(key=lambda x: datetime.strptime(re.search(r'_H(\d+)M(\d+)S(\d+)', x).group(0), '_H%HM%MS%S'))

    # Step 2: Find the earliest file for each index (0, 1, 2)
    earliest_files = {0: None, 1: None, 2: None}
    for file in Vk_files:
        index = int(file.split("_SingleShot(")[1].split(")_")[0])
        if earliest_files[index] is None:
            earliest_files[index] = file

    # Step 3: Group sorted files in sets of three, each containing one file with index 0, 1, 2
    Vk_files_grouped = []
    temp_group = {0: None, 1: None, 2: None}
    for file in Vk_files:
        index = int(file.split("_SingleShot(")[1].split(")_")[0])
        temp_group[index] = file
        if all(temp_group.values()):
            Vk_files_grouped.append([temp_group[0], temp_group[1], temp_group[2]])
            temp_group = {0: None, 1: None, 2: None}

    # Step 4: Calculate average g1 and standard deviation for each group of three files
    g1_values_per_group_with_delayiszero = []
    g1_std_per_group_with_delayiszero = []
    teff_values_per_group_with_delayiszero = []
    
    for files in Vk_files_grouped:
        if len(files) == 3:
            g1_values = []
            teff_values = []
            index = int(files[0].split("_SingleShot(")[1].split(")_")[0])
            ave_Vg_nc_file = earliest_files[index]  # Use the earliest file for the corresponding index for Vg_mean_I
            print(f"Using Vg file: {ave_Vg_nc_file}")
            Vg_mean_I = ave_Vg(ave_Vg_nc_file)

            for Vk_nc_file in files:
                print(f"Processing Vk file: {Vk_nc_file}")
                g1, teff = calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I, f01)
                g1_values.append(g1)
                teff_values.append(teff)

            # Average g1, standard deviation, and average teff for each group of three files
            g1_values_per_group_with_delayiszero.append(np.mean(g1_values))
            g1_std_per_group_with_delayiszero.append(np.std(g1_values))
            teff_values_per_group_with_delayiszero.append(np.mean(teff_values))
            
    g1_std_per_group=g1_std_per_group_with_delayiszero[1:]
    teff_values_per_group=teff_values_per_group_with_delayiszero[1:]
    # Step 5: Fit T1 using curve_fit from scipy for averaged g1 values
    
    g1_values_per_group = np.array(g1_values_per_group_with_delayiszero [1:])
    correlate_delays = np.array(correlate_delays[1:len(g1_values_per_group)+1])# Exclude the first data point
    A_guess = g1_values_per_group[0] - g1_values_per_group[-1]
    T1_guess = np.mean(correlate_delays)
    offset_guess = g1_values_per_group[-1]

    A_guess_upper = max(2 * (g1_values_per_group[0] - g1_values_per_group[-1]), 0.5 * (g1_values_per_group[0] - g1_values_per_group[-1]))
    A_guess_lower = min(2 * (g1_values_per_group[0] - g1_values_per_group[-1]), 0.5 * (g1_values_per_group[0] - g1_values_per_group[-1]))

    C_guess_upper = max(0.5 * g1_values_per_group[-1], 2 * g1_values_per_group[-1])
    C_guess_lower = min(0.5 * g1_values_per_group[-1], 2 * g1_values_per_group[-1])

    bounds = ((A_guess_lower, 0.1 * T1_guess, C_guess_lower), (A_guess_upper, 3 * T1_guess, C_guess_upper))
    p0 = (A_guess, T1_guess, offset_guess)

    ans, ans_error = curve_fit(t1, correlate_delays, g1_values_per_group, p0=p0, bounds=bounds)

    # Calculate g^(1) at delay = 0 using fitted model
    g1_at_zero_delay = t1(0, *ans)
    P_e = (1 / 2) - 1 / (2 * np.sqrt(1 + 4 * g1_at_zero_delay))
    qubit_frequency = 2 * np.pi * f01
    h_bar = 1.054571800e-34
    k_B = 1.381e-23
    effT = -h_bar * qubit_frequency / (k_B * np.log(P_e / (1 - P_e))) * 1000
    fitted_info = f"Fitted g^(1) at delay = 0: {g1_at_zero_delay:.4f}, Effective Temperature: {effT:.4f} mK"

    return correlate_delays, g1_values_per_group, g1_std_per_group, ans, fitted_info


def plot_all_results_in_one_figure(all_results):
    plt.figure(figsize=(12, 8))
    for plot_title, delays, g1_values, g1_std, ans, fitted_info in all_results:
        plt.errorbar(delays, g1_values, yerr=g1_std, marker='o', linestyle='none', capsize=5, label=f"{plot_title}")
        plt.plot(delays, t1(delays, *ans), linewidth=2)

    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("Average g^(1)")
    plt.legend()
    plt.grid(True)
    plt.title("Comparison of All Folders")
    plt.show()

def plot_each_result_individually(all_results):
    for plot_title, delays, g1_values, g1_std, ans, fitted_info in all_results:
        plt.figure(figsize=(10, 6))
        plt.errorbar(delays, g1_values, yerr=g1_std, marker='o', linestyle='none', capsize=5, label=f"{plot_title}")
        plt.plot(delays, t1(delays, *ans), linewidth=2, label="T1 Fit")
        plt.xlabel("Correlate Delay (us)")
        plt.ylabel("Average g^(1)")
        plt.legend()
        plt.grid(True)
        plt.title(plot_title)
        plt.show()

plot_option = input("Enter '1' to plot all results in one figure or '2' to plot each result individually: ")

all_results = []
for root, dirs, files in os.walk(root_directory):
    if 'SS' in dirs and 'SS_corre' in dirs:
        ss_folder = os.path.join(root, 'SS')
        ss_corre_folder = os.path.join(root, 'SS_corre')
        ss_files = [os.path.join(ss_folder, file) for file in os.listdir(ss_folder) if file.endswith('.nc')]
        if ss_files:
            ave_Ve_nc_file = ss_files[0]
            plot_title = os.path.basename(root)
            result = correlation_method(ave_Ve_nc_file, ss_corre_folder, f01, correlate_delays, plot_title)
            if result:
                all_results.append((plot_title, *result))

    # Plot all results in one figure
    plt.figure(figsize=(12, 8))
    for plot_title, delays, g1_values, g1_std, ans, fitted_info in all_results:
        plt.errorbar(delays, g1_values, yerr=g1_std, marker='o', linestyle='none', capsize=5, label=f"{plot_title}")
        plt.plot(delays, t1(delays, *ans), linewidth=2)

    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("Average g^(1)")
    plt.legend()
    plt.grid(True)
    plt.title("Comparison of All Folders")
    plt.show()

if __name__ == '__main__':
    root_directory = r"C:\Users\admin\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1026OrCouldBe1028"
    correlate_delays = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40]
    f01 = 4.31791e9
    
    if plot_option == '1':
        plot_all_results_in_one_figure(all_results)
    elif plot_option == '2':
        plot_each_result_individually(all_results)
    else:
        print("Invalid option. Please enter '1' or '2'.")
