import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from typing import Callable
from numpy import flip, pi, linspace, array, sqrt, std, mean, sort, diag, sign, absolute, arange
from Modularize.support import QDmanager, Data_manager
from qblox_drive_AS.support.Pulse_schedule_library import IQ_data_dis
from numpy import ndarray, cos, sin, deg2rad, real, imag, transpose, abs
from scipy.optimize import curve_fit

def plot_ROopti(Qmanager: QDmanager, nc_path: str,):
    # 检查 target_q 是否是列表
    ds = xr.open_dataset(nc_path)    
    for q in ds.data_vars:
        if q in Qmanager.refIQ:
            ref = Qmanager.refIQ[q]
        else:
            ref = [0, 0]
        # 进行绘图操作
        
        f, ro_amp = array(ds.coords['freq']), array(ds.coords['amp'])

        i_data, q_data = array(ds[q])[0].reshape(ro_amp.shape[0],f.shape[0]), array(ds[q])[1].reshape(ro_amp.shape[0],f.shape[0])
        
        # ro_amp = arange(ro_amp.shape[0])
        amp = sqrt((i_data - array(ref)[0])**2 + (q_data - array(ref)[1])**2).transpose()
        
        fig, ax1 = plt.subplots(figsize=(12, 8))
        c = ax1.pcolormesh(ro_amp**2, f * 1e-9, amp, cmap='viridis')#mV:(ro_amp*1e3)**2
        cbar =fig.colorbar(c, ax=ax1, pad=0.1)
        cbar.set_label(label='Contrast (V)',size=14)
        cbar.ax.tick_params(labelsize=14)
        
        for spine in ax1.spines.values():
            spine.set_linewidth(2)
        ax1.set_xlabel("Readout power ($V^2$)", fontsize=14)
        ax1.set_ylabel("Qubit Frequency (GHz)", fontsize=14)
        
        # 添加第二个 x 轴，共享 y 轴
        ax2 = ax1.twiny()
        # ax2.set_xlabel("Photon number", fontsize=14)
        # ax2.set_xlim(0, 5)  # 设置第二个 x 轴的范围

        # 添加第二个 y 轴，共享 x 轴
        ax3 = ax1.twinx()
        # ax3.set_ylabel("AC Stark Shift (MHz)", fontsize=14)
        # ax3.set_ylim(0, 40)  # 设置第二个 y 轴的范围

        ax1.xaxis.set_tick_params(labelsize=14,direction='in', width=2)
        ax1.yaxis.set_tick_params(labelsize=14,direction='in', width=2)
        ax2.xaxis.set_tick_params(labelsize=14,direction='in', width=2)
        ax3.yaxis.set_tick_params(labelsize=14,direction='in', width=2)
        
        fig.subplots_adjust(right=0.85)
        plt.title("AC Stark Shift", fontsize=20)
        plt.tight_layout()
        plt.show()
    

if __name__ == "__main__":
    QD_agent = QDmanager(r'C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_12_4\DRKE#242_SumInfo.pkl')
    QD_agent.QD_loader()
    file = r'C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_12_4\DRKEq0StarkshiftH14M20S7.nc'
    plot_ROopti(QD_agent, file)

   