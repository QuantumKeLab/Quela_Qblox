import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import mean, linspace
from numpy.random import randint

def t1(D, A, T1, offset):
    return A * np.exp(-D / T1) + offset

A=-1
C=0.5
t=np.linspace(0,50,100)
T=10

y = t1(t, A, T, C) + np.random.randint(0, 1, 100)

# Ensure lower bounds are less than upper bounds
A_guess_upper = max(2 * (y[0] - y[-1]), 0.5 * (y[0] - y[-1]))
A_guess_lower = min(2 * (y[0] - y[-1]), 0.5 * (y[0] - y[-1]))

C_guess_upper = max(0.5 * y[-1], 2 * y[-1])
C_guess_lower = min(0.5 * y[-1], 2 * y[-1])

ans, ans_error = curve_fit(t1, t, y, p0=(y[0] - y[-1], mean(t), y[-1]), bounds=((A_guess_lower, 0.1 * mean(t), C_guess_lower), (A_guess_upper, 3 * mean(t), C_guess_upper)))

plt.scatter(t, y, c='blue')
plt.plot(t, t1(t, *ans), c='red')
plt.grid()
plt.show()
