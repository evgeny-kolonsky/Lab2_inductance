# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:05:19 2023

@author: Gil
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.stats import linregress

#%%
xs = np.array([4,   4.2,   4.4,  4.6,  4.8,  5,    5.2,  5.4,  5.6,  5.8, 6, 6.2, 6.4, 6.6, 6.8, 7, 7.2, 7.4, 7.6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6, 9.8, 10, 10.5, 11, 11.5, 12, 12.5, 13])*1e-2 #m
Bs = np.array([130, 100.8, 83.1, 68.1, 55.3, 45.1, 38.4, 31.8, 27.8, 23.8, 21.0, 17.8, 15.9, 14.0, 12.5, 11.1, 9.96, 8.93, 7.71, 7.30, 6.57, 6.03, 5.50, 5.01, 4.64, 4.26, 3.95, 3.65, 3.41, 3.13, 2.93, 2.46, 2.11, 1.82, 1.60, 1.40,1.25  ])*1e-3 #T

xs -= 0.02 # correct origin to middle of magnet.

h= 2e-2 # m
b = 1e-2 #m

dist_dependence = (xs + h/2)/np.sqrt((xs + h/2)**2 + b**2) - (xs - h/2)/np.sqrt((xs - h/2)**2 + b**2)


reg = linregress(dist_dependence, Bs)
# print(reg)
M_reg = reg.slope
R_reg = reg.rvalue

inds = Bs<0.0075
reg2 = linregress(dist_dependence[inds], Bs[inds])
# print(reg)
M_reg2 = reg2.slope
R_reg2 = reg2.rvalue


plt.figure()
plt.plot(dist_dependence, Bs, '.', label = 'measurements')
plt.plot(dist_dependence, dist_dependence*reg.slope + reg.intercept, '-', label = 'reg full')
plt.plot(dist_dependence[inds], dist_dependence[inds]*reg2.slope + reg2.intercept, '-', label = 'reg B<0.0075T')
plt.grid()
plt.xlabel('distance [m]')
plt.ylabel('B [T]')
plt.ylim([0,0.0075])
plt.xlim([0,0.011])
plt.legend()



plt.figure()
plt.plot(xs, Bs, '.',label='measurements')
plt.plot(xs,dist_dependence*M_reg, '-', label='full range' )
plt.plot(xs,dist_dependence*M_reg2, '-', label='B < 0.0075T' )
plt.axhline(y=0.0075)
plt.grid()
plt.xlabel('distance [m]')
plt.ylabel('B [T]')

plt.legend()

#%%
