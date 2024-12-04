import matplotlib.pyplot as plt
import numpy as np
dis_mean_list=[3.4690, 4.2035,4.7803, 4.9041,5.0735, ]
dis_std_list=[0.03841,0.03792,0.04538,0.04781, 0.10572]

sigma_mean_list=[1.0110,0.9336,0.8573,0.8574,0.8513,]
sigma_std_list=[0.01028,0.00574,0.01057,0.00960,0.02420]
T_ro=np.linspace(1,2.00,len(dis_mean_list))
print(T_ro)

exp_snr=[10.7087,
13.2317,
14.571,
15.6795,
15.5059]
exp_snr_std=[0.06212,
0.1843,
0.30958,
0.38004,
0.41141]

snr=(np.array(dis_mean_list)/np.array(sigma_mean_list))*np.sqrt(T_ro)
snr_mean=np.log10(snr)*20
snr_std=snr_mean*np.sqrt((np.array(dis_std_list)/np.array(dis_mean_list))**2+(np.array(sigma_std_list)/np.array(sigma_mean_list))**2)
print('snr_mean=',snr_mean)
print('snr_std=',snr_std)
 
plt.figure(figsize=(8, 5))
plt.errorbar(T_ro, snr_mean, yerr=snr_std, fmt='o-', color='b', ecolor='r', capsize=5, label='SNR vs T_ro')
plt.errorbar(T_ro, exp_snr, yerr=exp_snr_std, fmt='o-', color='g', ecolor='purple', capsize=5, label='Experimental SNR vs T_ro')
 
plt.xlabel('T_ro (us)')
plt.ylabel('Power SNR (dB)')
plt.title('SNR as a function of T_ro with Error Bars')
plt.grid(True)
plt.legend()

 
plt.show()
