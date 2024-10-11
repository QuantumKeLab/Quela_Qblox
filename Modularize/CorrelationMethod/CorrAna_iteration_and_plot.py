import os
from xarray import open_dataset
from xarray import Dataset
import numpy as np
import matplotlib.pyplot as plt

def ave_Ve(ave_V_nc_file:str):
    dataset = open_dataset(ave_V_nc_file)
    e_var = dataset.variables['e'][:]
    e_I0 = e_var[0,:]
    Ve_mean_I = np.mean(e_I0)
    dataset.close()
    return Ve_mean_I 

def ave_Vg(ave_V_nc_file:str):
    dataset = open_dataset(ave_V_nc_file)
    g_var = dataset.variables['g'][:]
    g_I0 = g_var[0,:]
    Vg_mean_I = np.mean(g_I0)
    dataset.close()
    return Vg_mean_I 

def Vk_0_and_Vk_tau(Vk_nc_file:str):
    ds = open_dataset(Vk_nc_file)
    ss_dict = Dataset.to_dict(ds)
    pg_I = np.array(ss_dict['data_vars']['g']['data'][0])
    I_readout_1 = pg_I[0::2]
    I_readout_2 = pg_I[1::2]
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
    P_e = (1 / 2 )- 1 / (2 * np.sqrt(1 + 4 * normalized_g_1))
    qubit_frequency = 2 * np.pi * f01
    h_bar = 1.054571800e-34
    k_B = 1.381e-23
    T_eff = -h_bar * qubit_frequency / (k_B * np.log(P_e / (1 - P_e))) * 1000
    return normalized_g_1, T_eff

def main():
    # Step 1: Get Ve_mean_I
    ave_V_nc_file = r"C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_10\DRKEq1_cor_2\SS\DRKEq1_SingleShot(0)_H21M22S26.nc"
    Ve_mean_I = ave_Ve(ave_V_nc_file)
    Vg_mean_I = ave_Vg(ave_V_nc_file)
    # Step 2: List of Vk_nc_files and correlate_delay values
    Vk_directory = r"C:\Users\admin\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_10\DRKEq1_cor_2"
    Vk_nc_files = [os.path.join(Vk_directory, f) for f in os.listdir(Vk_directory) if f.endswith('.nc')]
    # Vk_nc_files.sort() # 確保按名稱排序
    # Vk_nc_files = Vk_nc_files[1:]
    correlate_delays = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    # correlate_delays = [0, 1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50]  # delays in µs
    # correlate_delays = list(range(len(Vk_nc_files)))  # Example: use file index as delay if no specific values are provided
    # correlate_delays = list(range(1, len(Vk_nc_files) + 1))

    f01=4.64e9
    print("qubit frequency: ", f01/1e9)
    # Step 3: Calculate g^(1) and T_eff for each Vk_nc_file
    g1_values = []
    teff_values = []
    for Vk_nc_file in Vk_nc_files:
        g1, teff = calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I,f01)
        g1_values.append(g1)
        teff_values.append(teff)

    # Step 4: Plot g^(1) vs correlate_delay
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(correlate_delays, g1_values, marker='o', color='b', label="g^(1)")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("g^(1)")
    # plt.title("g^(1) vs Correlate Delay")
    plt.legend()
    plt.grid(True)

    # Step 7: Plot Effective Temperature vs correlate_delay
    plt.subplot(2, 1, 2)
    plt.plot(correlate_delays, teff_values, marker='o', color='r', label="Effective Temperature (mK)")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("Effective Temperature (mK)")
    # plt.title("Effective Temperature vs Correlate Delay")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Optional: Print T_eff values for each delay
    # for delay, teff in zip(correlate_delays, teff_values):
    #     print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")

delay_taus = [1]*10
delay_taus = [2]*10

delay_taus = [5]*60

if __name__ == '__main__':
    main()
