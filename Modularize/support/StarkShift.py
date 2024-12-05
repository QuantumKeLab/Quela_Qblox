import os, json, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import xarray as xr
import quantify_core.data.handling as dh
import matplotlib.pyplot as plt
from typing import Callable
from numpy import flip, pi, linspace, array, sqrt, std, mean, sort, diag, sign, absolute, arange
from Modularize.support import QDmanager, Data_manager
from Modularize.support.Pulse_schedule_library import IQ_data_dis
from numpy import ndarray, cos, sin, deg2rad, real, imag, transpose, abs
from scipy.optimize import curve_fit
import numpy as np

# def plot_ROopti(Qmanager: QDmanager, nc_path: str,):
#     # 检查 target_q 是否是列表
#     ds = xr.open_dataset(nc_path)    
#     for q in ds.data_vars:
#         if q in Qmanager.refIQ:
#             ref = Qmanager.refIQ[q]
#         else:
#             ref = [0, 0]
#         # 进行绘图操作
        
#         f, ro_amp = array(ds.coords['freq']), array(ds.coords['amp'])

#         i_data, q_data = array(ds[q])[0].reshape(ro_amp.shape[0],f.shape[0]), array(ds[q])[1].reshape(ro_amp.shape[0],f.shape[0])
        
#         # ro_amp = arange(ro_amp.shape[0])
#         amp = sqrt((i_data - array(ref)[0])**2 + (q_data - array(ref)[1])**2).transpose()
        
#         fig, ax1 = plt.subplots(figsize=(12, 8))
#         c = ax1.pcolormesh(ro_amp**2, f * 1e-9, amp, cmap='viridis')#mV:(ro_amp*1e3)**2
#         cbar =fig.colorbar(c, ax=ax1, pad=0.1)
#         cbar.set_label(label='Contrast (V)',size=14)
#         cbar.ax.tick_params(labelsize=14)
        
#         for spine in ax1.spines.values():
#             spine.set_linewidth(2)
#         ax1.set_xlabel("Readout power ($V^2$)", fontsize=14)
#         ax1.set_ylabel("Qubit Frequency (GHz)", fontsize=14)
        
#         # 添加第二个 x 轴，共享 y 轴
#         ax2 = ax1.twiny()
#         # ax2.set_xlabel("Photon number", fontsize=14)
#         # ax2.set_xlim(0, 5)  # 设置第二个 x 轴的范围

#         # 添加第二个 y 轴，共享 x 轴
#         ax3 = ax1.twinx()
#         # ax3.set_ylabel("AC Stark Shift (MHz)", fontsize=14)
#         # ax3.set_ylim(0, 40)  # 设置第二个 y 轴的范围

#         ax1.xaxis.set_tick_params(labelsize=14,direction='in', width=2)
#         ax1.yaxis.set_tick_params(labelsize=14,direction='in', width=2)
#         ax2.xaxis.set_tick_params(labelsize=14,direction='in', width=2)
#         ax3.yaxis.set_tick_params(labelsize=14,direction='in', width=2)
        
#         fig.subplots_adjust(right=0.85)
#         plt.title("AC Stark Shift", fontsize=20)
#         plt.tight_layout()
#         plt.show()
    
def linear_fit(x, a, b):
    """线性函数用于拟合"""
    return a * x + b

def plot_ROopti(Qmanager, nc_path,dress_e,dress_g):
    ds = xr.open_dataset(nc_path)
    
    for q in ds.data_vars:
        if q in Qmanager.refIQ:
            ref = Qmanager.refIQ[q]
        else:
            ref = [0, 0]

        # 提取数据
        f, ro_amp = np.array(ds.coords['freq']), np.array(ds.coords['amp'])
        i_data = np.array(ds[q])[0].reshape(ro_amp.shape[0], f.shape[0])
        q_data = np.array(ds[q])[1].reshape(ro_amp.shape[0], f.shape[0])
        amp = np.sqrt((i_data - np.array(ref)[0])**2 + (q_data - np.array(ref)[1])**2).transpose()

        # 找出每列中最大值及其索引
        max_indices = np.argmax(amp, axis=0)  # 每列最大值的索引
        max_freqs = f[max_indices]  # 对应最大值的频率
        max_ro_power = ro_amp**2  # Readout Power    
        
        # 计算相对差值
        delta_freqs = max_freqs - max_freqs[0]

        # 打印 Readout Power = 0 和最大 Readout Power 的红色点数值
        print(f"Readout Power = 0: Max Frequency = {max_freqs[0]*1e-9:.3f} GHz")
        print(f"Max Readout Power = {max_ro_power[-1]:.3f} V^2: Max Frequency = {max_freqs[-1]*1e-9:.3f} GHz")


        """Parameter calculation"""
        AC_shift=(max_freqs[0]-max_freqs[-1])*1e-9
        print("AC Stark Shift (GHz):",AC_shift)
        
        # dress_g=5.85023
        # dress_e=5.84954
        # Kai_eff=dress_g-dress_e  #Need to get the info from ROcalibration
        # n_critical_shift=AC_shift/(2* Kai_eff) 
        # print("Critical photon number:",n_critical_shift)
        
        
    
        # dress_atte= 32 #Need to get the info from PowerDependence
        # V_dress=0.1#Need to get the info from PowerDependence
        # detuning=1#Need to get the info from Twotone
        # coupling_g=1#Need to get the info from Twotone
        
        # n_critical=(detuning/(2*coupling_g))**2
        # constant=10**(-dress_atte/10)*V_dress**2/n_critical
        # V_dress_new=(constant*n_critical/10**(-dress_atte/10))**(1/2)
        # print("New dress state readout voltage:", V_dress_new)
        
        
       
        
        """Fitting and plot"""
        # 线性拟合
        popt, pcov = curve_fit(linear_fit, max_ro_power, max_freqs)
        a, b = popt
        fit_line = linear_fit(max_ro_power, *popt)

        # 绘图
        fig, ax1 = plt.subplots(figsize=(12, 8))
        c = ax1.pcolormesh(max_ro_power, f * 1e-9, amp, cmap='viridis', shading='auto')
        cbar = fig.colorbar(c, ax=ax1, pad=0.1)
        cbar.set_label(label='Contrast (V)', size=14)
        cbar.ax.tick_params(labelsize=14)

        
        slope=abs(a)
        dress_g=5.85023
        dress_e=5.84954
        Kai_eff=dress_g-dress_e 
        A=slope/(2*Kai_eff*1e9)#= 8297.101
        max_ro_power=np.array(max_ro_power)

        n_photon=A*max_ro_power
        print(n_photon)
        
        # 绘制最大值点
        ax1.scatter(max_ro_power, max_freqs * 1e-9, color='red', label='Max Points', zorder=10)
        # 绘制拟合直线
        ax1.plot(max_ro_power, fit_line * 1e-9, 'k--', label=f'Linear Fit\ny = {a:.3e}x + {b:.3e}', linewidth=2)

        ax1.set_xlabel("Readout power ($V^2$)", fontsize=14)
        ax1.set_ylabel("Qubit Frequency (GHz)", fontsize=14)
        ax1.xaxis.set_tick_params(labelsize=14, direction='in', width=2)
        ax1.yaxis.set_tick_params(labelsize=14, direction='in', width=2)
        
        # ax2 = ax1.twinx()
        # ax2.set_ylabel("AC Stark Shift (Relative) [MHz]", fontsize=14)
        # ax2.set_ylim(delta_freqs[-1] * 1e-6, 0)  # 转为 MHz
        # ax2.yaxis.set_tick_params(labelsize=14, direction='in', width=2)

        ax3 = ax1.twiny()
        ax3.set_xlabel("Photon number", fontsize=14)
        ax3.set_xlim(0, np.max(n_photon))  # 使用 n_photon 的最大值
        ax3.xaxis.set_tick_params(labelsize=14, direction='in', width=2)

        
        # 添加图例
        ax1.legend(fontsize=12, loc='upper right')
        
    

        # 添加标题和调整
        plt.title("AC Stark Shift with Linear Fit", fontsize=20)
        plt.tight_layout()
        plt.show()

        # 打印拟合参数
        print(f"Linear fit parameters: a = {a:.3e}, b = {b:.3e}")
        print(f"Covariance matrix:\n{pcov}")
        
         # 打印拟合参数
        print(f"Linear fit parameters: a = {a:.3e}, b = {b:.3e}")
        print(f"Covariance matrix:\n{pcov}")

        """绘制微分图"""
        # 计算微分
        d_freqs = np.gradient(max_freqs, max_ro_power)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(max_ro_power, d_freqs * 1e-9, 'b-o', label='Frequency Gradient (GHz/$V^2$)')
        ax.set_xlabel("Readout power ($V^2$)", fontsize=14)
        ax.set_ylabel("Gradient (GHz/$V^2$)", fontsize=14)
        ax.xaxis.set_tick_params(labelsize=12, direction='in', width=2)
        ax.yaxis.set_tick_params(labelsize=12, direction='in', width=2)
        ax.legend(fontsize=12)
        plt.title("Frequency Gradient vs Readout Power", fontsize=16)
        plt.tight_layout()
        plt.show()
        
        return AC_shift#, n_critical_shift, V_dress_new

if __name__ == "__main__":
    QD_agent = QDmanager(r'C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\QD_backup\2024_12_4\DRKE#242_SumInfo.pkl')
    QD_agent.QD_loader()
    file = r"C:\Users\Ke Lab\Documents\GitHub\Quela_Qblox\Modularize\Meas_raw\2024_12_5\DRKEq0StarkshiftH19M23S42.nc"
    dress_g=5.85023
    dress_e=5.84954
    plot_ROopti(QD_agent, file,dress_e,dress_g)
    

   