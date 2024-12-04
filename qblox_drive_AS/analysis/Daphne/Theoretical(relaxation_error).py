import numpy as np
import matplotlib.pyplot as plt

# 定義參數
T1 = 8.8 # 設定 T1 值
T_ro = np.linspace(0.75, 2, 6)  # 時間範圍 np.linspace(0, 5 * T1, 500)
p10_exp_list=[11.5071,
12.6858,
15.7801,
16.5618,
18.4811,
19.4587]

p10_exp_std_list=[2.3447,
1.19023,
0.68147,
0.46149,
0.81586,
1.8515]
# 計算函數值
F_T = (1/T_ro) * np.array([np.trapz(1 - np.exp(-np.linspace(0, t, 1000) / T1), np.linspace(0, t, 1000)) if t > 0 else 0 for t in T_ro])

# 繪圖
plt.figure(figsize=(8, 5))
plt.plot(T_ro, F_T, label=r"$\frac{1}{T} \int_0^T (1 - e^{-t/T_1}) dt$")
plt.errorbar(T_ro, p10_exp_list, yerr=p10_exp_std_list, fmt='o-', color='g', ecolor='purple', capsize=5, label='Experimental SNR vs T_ro')
# plt.axhline(1, color="red", linestyle="--", label="Asymptote at 1")
plt.xlabel(r"$T$")
plt.ylabel("Value")
# plt.title(r"Plot of $\frac{1}{T} \int_0^T (1 - e^{-t/T_1}) dt$")
plt.legend()
plt.grid()
plt.show()

""""""

# # 定義參數範圍
# T1_values = [8.0, 8.8, 9.6]  # 下界、中值、上界
# colors = ['green', 'blue', 'orange']  # 曲線顏色
# labels = [r"$T_1 = 8.0$", r"$T_1 = 8.8$", r"$T_1 = 9.6$"]

# # 時間範圍
# T = np.linspace(0, 5 * max(T1_values), 500)

# # 繪圖
# plt.figure(figsize=(8, 5))
# for T1, color, label in zip(T1_values, colors, labels):
#     # 計算函數值
#     F_T = (1 / T) * np.array([
#         np.trapz(1 - np.exp(-np.linspace(0, t, 1000) / T1), np.linspace(0, t, 1000)) if t > 0 else 0 for t in T
#     ])
#     # 繪製曲線
#     plt.plot(T, F_T, label=label, color=color)

# # 添加極限虛線
# plt.axhline(1, color="red", linestyle="--", label="Asymptote at 1")

# # 添加標籤和標題
# plt.xlabel(r"$T$ (time)")
# plt.ylabel("Value")
# plt.title(r"Plot of $\frac{1}{T} \int_0^T (1 - e^{-t/T_1}) dt$ with $T_1 = 8.8 \pm 0.8$")
# plt.legend()
# plt.grid()
# plt.show()

