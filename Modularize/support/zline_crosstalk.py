import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..\..'))
import matplotlib.pyplot as plt
from numpy import sqrt
from Modularize.support.QDmanager import QDmanager
from Modularize.support.QuFluxFit import convert_netCDF_2_arrays
from qcat.zline_crosstalk.ramsey_2dfft import analysis_crosstalk_value
from qcat.visualization.zline_crosstalk import plot_analysis

QD_agent = QDmanager(r"D:\SynologyDrive\12 Codes\RF\20240224_QM_opt\Quela_Qblox\Modularize\QD_backup\2024_5_10\DRKE#116_SumInfo.pkl")
QD_agent.QD_loader()
nc_path = r"d:\SynologyDrive\09 Data\Fridge Data\Qubit\20240501_DRKe_5Q4C#6\20240509\q4\Zline_Crosstalk\q4z3\DRKEq0_Flux2tone_H11M3S55.nc"
ref = QD_agent.refIQ['q1']

d_z_target_amp, d_z_crosstalk_amp, i, q = convert_netCDF_2_arrays(nc_path)
data = sqrt((i - i[len(i)//2][len(i)//2])**2 + (q - q[len(q)//2][len(q)//2])**2).transpose()

analysis_crosstalk_value(d_z_crosstalk_amp, d_z_target_amp, data.transpose())
plot_analysis(d_z_crosstalk_amp, d_z_target_amp, data.transpose())