import xarray as xr
import quantify_core.data.handling as dh
from Modularize.support.QDmanager import QDmanager
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from numpy import sqrt, array, cos, meshgrid, pi
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr
from Modularize.support import init_meas


# from quantify_core.analysis.spectroscopy_analysis import ResonatorSpectroscopyAnalysis
# from quantify_core.analysis.base_analysis import Basic2DAnalysis



def plot_QbFlux(Qmanager:QDmanager, nc_path:str, target_q:str):
    ref = Qmanager.refIQ[target_q]
    # plot flux-qubit 
    f,z,i,q = convert_netCDF_2_arrays(nc_path)
    amp = array(sqrt((i-array(ref)[0])**2+(q-array(ref)[1])**2)).transpose()
    fig, ax = plt.subplots()
    c = ax.pcolormesh(z, f, amp, cmap='RdBu')
    fig.colorbar(c, ax=ax)
    plt.show()


def plot_Zline_Crosstalk(Qmanager: QDmanager, nc_path: str, target_q: str):
    ref = Qmanager.refIQ[target_q]
    # plot flux-qubit 
    f, z, i, q = convert_netCDF_2_arrays(nc_path)
    amp = sqrt((i - array(ref)[0])**2 + (q - array(ref)[1])**2).transpose()

    # 定義要拟合的函數
    def func(data, A, B, C, D, E):
        z, f = data
        return A * cos(B * z + C * f + D) + E

    # 整理數據
    z_flat = z.flatten()
    f_flat = f.flatten()
    amp_flat = amp.flatten()

    # 生成 z 和 f 的所有配對組合
    Z_pairs, F_pairs = meshgrid(z_flat, f_flat)

    # 使用 curve_fit 進行拟合
    popt, pcov = curve_fit(func, (Z_pairs.flatten(), F_pairs.flatten()), amp_flat)

    # 輸出拟合參數
    print("拟合参数:", popt)

    # 繪製圖像
    fig, ax = plt.subplots()
    c = ax.pcolormesh(z, f, amp, cmap='RdBu')
    fig.colorbar(c, ax=ax)

    # 繪製拟合曲線
    Z, F = meshgrid(z_flat, f_flat)
    B, C, D = popt[1:-1]  # 不包括 A, E

    # 計算 B*z + C*f + D = 0 和 B*z + C*f + D = ±π/2
    fit_zero = (-D - B * Z) / C
    fit_plus_pi_2 = (pi / 2 - D - B * Z) / C
    fit_minus_pi_2 = (-pi / 2 - D - B * Z) / C

    # 繪製等值線
    ax.contour(Z, F, fit_zero, levels=[0], colors='k')
    ax.contour(Z, F, fit_plus_pi_2, levels=[0], colors='b', linestyles='--')
    ax.contour(Z, F, fit_minus_pi_2, levels=[0], colors='b', linestyles='--')

    plt.show()




# from quantify_scheduler.helpers.collections import find_port_clock_path
# QD_agent = QDmanager('Modularize/QD_backup/2024_3_19/DR2#171_SumInfo.pkl')
# QD_agent.QD_loader()
# print(QD_agent.Fluxmanager.get_sweetBiasFor('q2'))
# print(QD_agent.Fluxmanager.get_PeriodFor('q2'))
# qd = QD_agent.quantum_device
# hcfg = QD_agent.Hcfg
# qubit = qd.get_element('q4')
# output_path = find_port_clock_path(
#         hcfg, port=qubit.ports.readout(), clock=qubit.name + ".ro"
#     )

# cluster_key, module_key, output_key, _, _ = tuple(output_path)
# readout_module = hcfg[cluster_key][module_key]
# print(readout_module[output_key]["output_att"])


# meas_datadir = 'tempt'
# dh.set_datadir(meas_datadir)
# print(ds.attrs["tuid"])
# ana = Basic2DAnalysis(tuid=ds.attrs["tuid"], dataset=ds).run()
# print(ana.quantities_of_interest)



# from Modularize.support import QDmanager
# from numpy import array
# QD_path = 'Modularize/QD_backup/2024_3_31/DR2#171_SumInfo.pkl'
# QD_agent = QDmanager(QD_path)
# QD_agent.QD_loader()
# for i in ["q0","q1","q2","q3","q4"]:
#     qu = QD_agent.quantum_device.get_element(i)
#     rof = qu.clock_freqs.readout()
#     xyf = qu.clock_freqs.f01()
#     xyl = qu.rxy.amp180()
#     T1 = QD_agent.Notewriter.get_T1For(target_q=i)
#     T2 = QD_agent.Notewriter.get_T2For(target_q=i)
#     print(f"***{i}:")
#     print(f"rof : {rof}")
#     print(f"xyf : {xyf}")
#     print(f"xyl : {xyl}")
#     print(f"T1 : {T1}")
#     print(f"T2 : {T2}")
#     print("===========================\n")
# def aprx_fq(disper_MHz:float,bareF_MHz:float,g_MHz:float=45.0):
#     return bareF_MHz-(g_MHz**2)/disper_MHz


# print("aprx fq=",aprx_fq(x,bare))

if __name__ == '__main__':
    # QD_agent = QDmanager('Modularize/QD_backup/2024_5_2/DR1#11_SumInfo.pkl')
    # QD_agent.QD_loader()
    # print(QD_agent.quantum_device.get_element("q0").clock_freqs.f01())
    DRandIP = {"dr":"drke","last_ip":"116"}
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    QD_agent, cluster, meas_ctrl, ic, Fctrl = init_meas(QuantumDevice_path=QD_path,mode='l')
    qb=''
    qubit = QD_agent.quantum_device.get_element(qb)
    qubit.clock_freqs.readout(FD_results[qb].quantities_of_interest["freq_0"])
    QD_agent.Fluxmanager.save_sweetspotBias_for(target_q=qb,bias=FD_results[qb].quantities_of_interest["offset_0"].nominal_value)
