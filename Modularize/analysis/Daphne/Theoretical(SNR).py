import matplotlib.pyplot as plt
import numpy as np
dis_mean_list=[3.4690, 4.2035,4.7803, 4.9041,5.0735, ]
dis_std_list=[0.03841,0.03792,0.04538,0.04781, 0.10572]

sigma_mean_list=[1.0110,0.9336,0.8573,0.8574,0.8513,]
sigma_std_list=[0.01028,0.00574,0.01057,0.00960,0.02420]
T_ro=np.linspace(1,2.00,len(dis_mean_list))
print(T_ro)

snr_mean=(np.array(dis_mean_list)/np.array(sigma_mean_list))*np.sqrt(T_ro)
snr_std=snr_mean*np.sqrt((np.array(dis_std_list)/np.array(dis_mean_list))**2+(np.array(sigma_std_list)/np.array(sigma_mean_list))**2)
print('snr_mean=',snr_mean)
 
plt.figure(figsize=(8, 5))
plt.errorbar(T_ro, snr_mean, yerr=snr_std, fmt='o-', color='b', ecolor='r', capsize=5, label='SNR vs T_ro')

 
plt.xlabel('T_ro')
plt.ylabel('SNR')
plt.title('SNR as a function of T_ro with Error Bars')
plt.grid(True)
plt.legend()

 
plt.show()
