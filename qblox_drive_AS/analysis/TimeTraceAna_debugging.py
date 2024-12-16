import os, sys, json 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from xarray import Dataset, open_dataset
from qblox_drive_AS.support import QDmanager, Data_manager
import matplotlib.pyplot as plt
from numpy import ndarray, array, median, std, mean
from qblox_drive_AS.support.Pulse_schedule_library import IQ_data_dis, dataset_to_array, T2_fit_analysis, T1_fit_analysis
from datetime import datetime 
from matplotlib.gridspec import GridSpec as GS
from Modularize.analysis.Radiator.RadiatorSetAna import sort_set
import numpy as np

def time_label_sort(nc_file_name: str):
    return datetime.strptime(nc_file_name.split("_")[-1].split(".")[0], "H%HM%MS%S")

def plot_coherence_timetrace(raw_data: ndarray, time_samples: ndarray, ans: ndarray, q: str, raw_data_folder: str, exp: str, detunings: ndarray = []):
    """
    Plot color map, histogram, and detuning if ramsey in one figure
    """
    json_files = [os.path.join(raw_data_folder, name) for name in os.listdir(raw_data_folder) if (os.path.isfile(os.path.join(raw_data_folder, name)) and name.split(".")[-1] == "json")]

    if json_files:
        time_json_path = json_files[0]
        with open(time_json_path) as time_record_file:
            time_past_dict: dict = json.load(time_record_file)
        time_array = array(list(time_past_dict.values())[0])
    else:
        print("No JSON files found in the specified folder.")
        return

    median_ans, std_ans = median(ans), std(ans)
    fig = plt.figure(dpi=100, figsize=(12, 9))
    
    gs = GS(2, 2, width_ratios=[2, 1], height_ratios=[1, 5])
    if exp.lower() == 't2':
        ax0 = fig.add_subplot(gs[0, 0])
        d = ax0.plot(time_array, detunings, c='magenta')
        ax0.set_xlim(min(time_array), max(time_array))
        ax0.set_title("Transition detuning (MHz)", fontsize=20)
        ax0.grid()
        ax0.xaxis.set_tick_params(labelsize=16)
        ax0.yaxis.set_tick_params(labelsize=16)

    ax1 = fig.add_subplot(gs[:, 1])
    ax1.hist(ans, bins='auto', density=False)
    ax1.axvline(median_ans, c='k', ls='--', lw='1')
    ax1.set_ylabel("Counts", fontsize=20)
    ax1.set_xlabel(f"{exp.upper()} (us)", fontsize=20)
    ax1.set_title(f"{exp.upper()} = {round(median_ans, 1)} $\pm$ {round(std_ans, 1)} us", fontsize=20)

    if exp.lower() in ['t1', 'os']:
        ax2 = fig.add_subplot(gs[:, 0])
    elif exp.lower() == 't2':
        ax2 = fig.add_subplot(gs[1, 0])
    else:
        pass
    c = ax2.pcolormesh(time_array, time_samples * 1e6, raw_data.transpose() * 1000, cmap='RdBu')
    fig.colorbar(c, ax=ax2, label='Contrast (mV)', location='bottom', )
    ax2.plot(time_array, ans, c='green')
    ax2.axhline(median_ans + std_ans, linestyle='--', c='orange')
    ax2.axhline(median_ans, linestyle='-', c='orange')
    ax2.axhline(median_ans - std_ans, linestyle='--', c='orange')
    ax2.set_xlabel('Time past (min)', fontsize=20)
    ax2.set_ylabel("Free evolution time (us)", fontsize=20)

    for Ax in [ax1, ax2]:
        Ax: plt.Axes
        Ax.xaxis.set_tick_params(labelsize=16)
        Ax.yaxis.set_tick_params(labelsize=16)

    plt.title(f"Time dependent {exp.upper()}", fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(raw_data_folder, f"{exp.upper()}_{q}_timeDep.png"))
    plt.close()
    
def colormap(x:ndarray, y:ndarray, z:ndarray, fit_values:ndarray, ax:plt.Axes=None, fig_path:str=None):
    if ax is None:
        fig, ax = plt.subplots()
    c = ax.pcolormesh(x, y, z.transpose(), cmap='RdBu')
    ax.plot(x,fit_values, c='green')
    ax.axhline(median(fit_values)+std(fit_values),linestyle='--',c='orange')
    ax.axhline(median(fit_values),linestyle='-',c='orange')
    ax.axhline(median(fit_values)-std(fit_values),linestyle='--',c='orange')
    ax.set_xlabel('Time past (min)',fontsize=20)
    ax.set_ylabel("Free evolution time (us)",fontsize=20)
    plt.colorbar(c, ax=ax, label='Contrast (V)')
    plt.tight_layout()
    if fig_path is not None:
        plt.savefig(fig_path)
    plt.close()

def plot_timeDepCohe(time_values:ndarray, y_values:ndarray, exp:str, fig_path:str=None, units:dict={"x":"min","y":"µs"}):
    fig, axs = plt.subplots(1,3,figsize=(16,8))
    ax:plt.Axes = axs[0]
    ax.grid()
    ax.plot(time_values, y_values)
    ax.axhline(median(y_values)+std(y_values),linestyle='--',c='orange',label='1σ')
    ax.axhline(mean(y_values),linestyle='-',c='orange',label='mean')
    ax.axhline(median(y_values)-std(y_values),linestyle='--',c='orange')
    ax.axhline(median(y_values),linestyle='-.',c='red',label='median')
    ax.set_xlabel(f"time past ({units['x']})",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})",fontsize=16)
    ax.set_title(f"Time dependent {exp.upper()}",fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.legend()

    ax:plt.Axes=axs[1]
    ax.grid()
    ax.hist(y_values, bins='auto', density=False,orientation="horizontal")
    ax.axhline(median(y_values),c='k',ls='--',lw='1')
    ax.set_xlabel("Counts",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})", fontsize=16)
    ax.set_title(f"{exp.upper()} = {round(median(y_values),2 if exp == 'δf' else 1)} $\pm$ {round(std(y_values),2 if exp == 'δf' else 1)} {units['y']}",fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)

    ax:plt.Axes=axs[2]
    ax.hist(y_values, bins=5000, orientation="horizontal", cumulative=True, density=True,histtype="step")
    ax.axhline(median(y_values),c='k',ls='--',lw='1')
    ax.axvline(0.5,c='red',ls='--',lw='1',label='density=0.5')
    ax.set_xlabel("Density",fontsize=16)
    ax.set_ylabel(f"{exp.upper()} ({units['y']})", fontsize=16)
    ax.set_title("CDF",fontsize=20)
    ax.legend()
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)

    plt.tight_layout(pad=5.0)
    if fig_path is not None:
        plt.savefig(fig_path)
    plt.close()


# 其他函數保持不變

if __name__ == "__main__":
    folder_paths = {"T1_folder_path": r"C:\Users\Ke Lab\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\2024_10_27\T1_overnight"}
    QD_file_path = r"C:\Users\Ke Lab\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\QD_backup\2024_10_26\DRKE#242_SumInfo.pkl"
    qs = ['q0']
    qs = 'q0'
    sort_mode = 'time'
    T1_guess_parameter = 10e-6
    T2_guess_parameter = 15e-6

    QD_agent = QDmanager(QD_file_path)
    QD_agent.QD_loader()

    for folder_name in folder_paths:
        folder = folder_paths[folder_name]
        if folder_paths[folder_name] != '':
            if sort_mode == 'time':
                files = sorted([name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder, name)) and name.split(".")[-1] == "nc")], key=lambda name: time_label_sort(name))
            elif sort_mode == 'idx':
                files = [name for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder, name)) and name.split(".")[-1] == "nc")]
                sort_set(files, 3)
            else:
                raise KeyError(f"Unsupported sort mode was given = '{sort_mode}'")

            raw_data = []
            ans = []
            detu = []

            for idx, file in enumerate(files):
                path = os.path.join(folder, file)
                nc = open_dataset(path)
                if folder_name.split("_")[0] in ["T1", "T2"]:
                    samples = array(nc.x0)
                    I, Q = dataset_to_array(dataset=nc, dims=1)
                    data = IQ_data_dis(I, Q, ref_I=QD_agent.refIQ[qs][0], ref_Q=QD_agent.refIQ[qs][-1])
                    raw_data.append(data)

                    if folder_name.split("_")[0] == "T1":
                        data_fit = T1_fit_analysis(data=data, freeDu=samples, T1_guess=T1_guess_parameter)
                        print(data_fit.attrs['T1_fit'])
                        t1_fit = data_fit.attrs['T1_fit']
                        if not np.isnan(t1_fit) and 0 < t1_fit < 50:
                            ans.append(t1_fit * 1e6)
                        else:
                            if ans:
                                mean_ans = np.mean(ans)
                                std_ans = np.std(ans)
                                ans.append(mean_ans + std_ans)
                    else:
                        data_fit = T2_fit_analysis(data=data, freeDu=samples, T2_guess=T2_guess_parameter)
                        t2_fit = data_fit.attrs['T2_fit']
                        if not np.isnan(t2_fit):
                            ans.append(t2_fit * 1e6)
                        else:
                            if ans:
                                mean_ans = np.mean(ans)
                                std_ans = np.std(ans)
                                ans.append(mean_ans + std_ans)

                        if not np.isnan(data_fit.attrs['f']):
                            detu.append(data_fit.attrs['f'] * 1e-6)
                else:
                    pass
            print(f"Number of processed files: {len(files)}")

            json_files = [os.path.join(folder, name) for name in os.listdir(folder) if (os.path.isfile(os.path.join(folder, name)) and name.split(".")[-1] == "json")]

            if json_files:
                time_json_path = json_files[0]
                with open(time_json_path) as time_record_file:
                    time_past_dict = json.load(time_record_file)
                time_array = np.array(list(time_past_dict.values())[0])
                print(f"Length of time array: {len(time_array)}")

                if folder_name.split("_")[0].lower() in ['t1', 't2']:
                    colormap(time_array, samples * 1e6, np.array(raw_data), np.array(ans), fig_path=os.path.join(folder, f"{qs}_{folder_name.split('_')[0]}_timeDep_colormap.png"))

                if folder_name.split("_")[0].lower() == 't1':
                    plot_timeDepCohe(time_array, np.array(ans), folder_name.split("_")[0], units={"x": "min", "y": "µs"}, fig_path=os.path.join(folder, f"{qs}_T1_timeDep.png"))
                elif folder_name.split("_")[0].lower() == 't2':
                    plot_timeDepCohe(time_array, np.array(ans), folder_name.split("_")[0], units={"x": "min", "y": "µs"}, fig_path=os.path.join(folder, f"{qs}_T2_timeDep.png"))
                    plot_timeDepCohe(time_array, np.array(detu) - np.array(detu)[0], "δf", units={"x": "min", "y": "MHz"}, fig_path=os.path.join(folder, f"{qs}_Detune_timeDep.png"))
                else:
                    pass
            else:
                print("No JSON files found in the specified folder.")