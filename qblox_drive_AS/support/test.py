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
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from numpy import array, sqrt

# 假设有示例数据
f = np.linspace(0, 10, 100)  # frequency
ro_amp = np.linspace(0, 20, 50)  # readout amplitude
X, Y = np.meshgrid(ro_amp, f)
Z = np.sin(X) * np.cos(Y)  # 示例数据，用于展示的 color map

# 绘制 colormap
fig, ax1 = plt.subplots(figsize=(12, 8))
c = ax1.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
cbar = fig.colorbar(c, ax=ax1, label='Contrast (V)', pad=0.1)

ax1.set_xlabel("Readout power ($V^2$)", fontsize=16)
ax1.set_ylabel("Qubit Frequency (GHz)", fontsize=16)

# 找到每一列中的最大值的位置
max_indices = np.argmax(Z, axis=0)  # 找到每一列的最大值索引

# 获取最大值的坐标
max_ro_amp = ro_amp  # x 轴坐标是 readout amplitude 的值
max_frequencies = f[max_indices]  # y 轴坐标是最大值对应的频率

# 在 color map 上标出每一个最大值点
ax1.scatter(max_ro_amp, max_frequencies, color='red', s=50, edgecolor='black', label='Max Points')

# 可选：给每个点添加标注
for i in range(len(max_ro_amp)):
    ax1.annotate(f"({max_ro_amp[i]:.2f}, {max_frequencies[i]:.2f})",
                 (max_ro_amp[i], max_frequencies[i]),
                 textcoords="offset points",
                 xytext=(0, 5),  # 标注文本相对于点的偏移
                 ha='center',
                 fontsize=8)

# 添加图例
ax1.legend()

# 调整边距和显示
fig.subplots_adjust(right=0.85)
plt.title("AC Stark Shift with Max Points", fontsize=20)
plt.tight_layout()
plt.show()
