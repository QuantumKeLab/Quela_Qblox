import sys, os 
sys.path.append(r"C:\Users\User\Documents\GitHub\Quela_Qblox")
import xarray as xr
from xarray import open_dataset
from xarray import Dataset
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
from Modularize.support.Pulse_schedule_library import T1_fit_analysis

def ave_Ve(ave_V_nc_file:str):
    dataset = open_dataset(ave_V_nc_file)
    e_var = dataset.variables['e'][:]
    e_I0 = e_var[0,:]
    Ve_mean_I = np.mean(e_I0)
    dataset.close()
    return Ve_mean_I 

# def ave_Vg(ave_V_nc_file:str):
#     dataset = open_dataset(ave_V_nc_file)
#     g_var = dataset.variables['g'][:]
#     g_I0 = g_var[0,:]
#     Vg_mean_I = np.mean(g_I0)
#     dataset.close()
#     return Vg_mean_I 

def ave_Vg(ave_Vg_nc_file:str):
    ds = open_dataset(ave_Vg_nc_file)
    ss_dict = Dataset.to_dict(ds)
    pg_I = np.array(ss_dict['data_vars']['g']['data'][0])
    I_readout_1 = pg_I[0::2]
    Vg_mean_I = np.mean(I_readout_1)
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
    print(f"sum_result{sum_result}")
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

def correlation_method(ave_Vg_nc_file,ave_Ve_nc_file,Vk_directory,f01,correlate_delays): 
    # Step 1: Get Ve_mean_I
    Ve_mean_I = ave_Ve(ave_Ve_nc_file)
    Vg_mean_I = ave_Vg(ave_Vg_nc_file)
    # Step 2: List of Vk_nc_files and correlate_delay values
    Vk_nc_files = [os.path.join(Vk_directory, f) for f in os.listdir(Vk_directory) if f.endswith('.nc')]
    # Vk_nc_files.sort() # 確保按名稱排序
    # Vk_nc_files = Vk_nc_files[1:]
    # correlate_delays = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    # correlate_delays = [0, 1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50]  # delays in µs
    # correlate_delays = list(range(len(Vk_nc_files)))  # Example: use file index as delay if no specific values are provided
    # correlate_delays = list(range(1, len(Vk_nc_files) + 1))
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

#     # Step 7: Plot Effective Temperature vs correlate_delay
#     plt.subplot(2, 1, 2)
#     plt.plot(correlate_delays, teff_values, marker='o', color='r', label="Effective Temperature (mK)")
#     plt.xlabel("Correlate Delay (us)")
#     plt.ylabel("Effective Temperature (mK)")
#     # plt.title("Effective Temperature vs Correlate Delay")
#     plt.legend()
#     plt.grid(True)
#     plt.show()

    #Optional: Print T_eff values for each delay
    for delay, teff in zip(correlate_delays, teff_values):
        print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")
    return

# Assuming T1_func and T1_func_model are defined elsewhere
# This is how we would typically define an exponential function for fitting
def T1_func(D, A, T1, offset):
    return A * np.exp(-D / T1) + offset

# Assuming this is your model for fitting
T1_func_model = Model(T1_func)

# Example function to analyze correlation data
def correlation_method_fitting(ave_Ve_nc_file, ave_Vg_nc_file, Vk_directory, f01, correlate_delays): 
    # Step 1: Get Ve_mean_I
    Ve_mean_I = ave_Ve(ave_Ve_nc_file)
    Vg_mean_I = ave_Vg(ave_Vg_nc_file)
    
    # Step 2: List of Vk_nc_files and correlate_delay values
    Vk_nc_files = [os.path.join(Vk_directory, f) for f in os.listdir(Vk_directory) if f.endswith('.nc')]

    # Step 3: Calculate g^(1) and T_eff for each Vk_nc_file
    g1_values = []
    teff_values = []
    for Vk_nc_file in Vk_nc_files:
        g1, teff = calculate_g1_and_teff(Vk_nc_file, Vg_mean_I, Ve_mean_I, f01)
        g1_values.append(g1)
        teff_values.append(teff)

    # Convert g1_values and correlate_delays to np arrays
    g1_values = np.array(g1_values)
    correlate_delays = np.array(correlate_delays)

    # # Step 4: Plot g^(1) vs correlate_delay
    # plt.figure(figsize=(10, 6))
    # plt.subplot(2, 1, 1)
    # plt.plot(correlate_delays, g1_values, marker='o', color='b', label="g^(1)")
    # plt.xlabel("Correlate Delay (us)")
    # plt.ylabel("g^(1)")
    # plt.legend()
    # plt.grid(True)
    
    # Step 5: Perform the T1 exponential fitting
    T1_guess = 10e-6  # You can change this guess based on prior knowledge
    result_ds = T1_fit_analysis(g1_values, correlate_delays, T1_guess=T1_guess)

    # Step 6: Plot the fitted curve
    plt.subplot(2, 1, 2)
    plt.plot(result_ds['freeDu'], g1_values, 'bo', label="Data")
    plt.plot(result_ds['para_fit'], result_ds['fitting'], 'r-', label="Fitted curve")
    plt.xlabel("Correlate Delay (us)")
    plt.ylabel("g^(1)")
    plt.legend()
    plt.grid(True)

    # Step 7: Show the plots
    plt.tight_layout()
    plt.show()
     #Optional: Print T_eff values for each delay
    for delay, teff in zip(correlate_delays, teff_values):
        print(f"Delay {delay} µs: Effective Temperature = {teff:.2f} mK")
    # Return the analysis result
    return result_ds

# Example usage:
# analysis_result = correlation_method_fitting('ave_Ve.nc', 'ave_Vg.nc', '/path/to/Vk_files/', f01, correlate_delays)

if __name__ == '__main__':
    ave_Ve_nc_file = r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_24\13mK\SS\DRKEq0_SingleShot(0)_H18M47S9.nc"
    ave_Vg_nc_file = r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_24\13mK\SS_corre\DRKEq0_SingleShot(0)_H19M1S8.nc"
    
    Vk_directory = r"C:\Users\User\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_10_24\13mK\SS_corre"
    # correlate_delays = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30] 
    # correlate_delays = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    correlate_delays = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,4,6,8,10,15,20,25,30,35,40]   # delays in µs
    f01=4.31791e9
    correlation_method_fitting(ave_Ve_nc_file,ave_Vg_nc_file, Vk_directory,f01,correlate_delays)
