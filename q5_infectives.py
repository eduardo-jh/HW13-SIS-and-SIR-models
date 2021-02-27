#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW13 - Questions 5 and 6. SIS and SIR models

Created on Sat Feb 27 03:03:25 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

# d = 200  # time of infection, hours
d = 20  # question 6
alpha = 0.001  #rate of infection 1/person-hr
N = 400  # total population
B = 1/d  # rate of recovery
steps = 30  # time of simulation
dt = 0.1  # time step, hours
K = (N - B/alpha)
r = K*alpha
I0 = 1
S0 = N - I0

t = np.linspace(0, steps, int(steps/dt)+1)

I = np.zeros(len(t))
S = np.zeros(len(t))
I[0], S[0] = I0, S0

# Numerical solution, using Euler method
for i in range(1, len(t)):
    I[i] = I[i-1] + (alpha * I[i-1]*S[i-1] - B*I[i-1]) * dt
    S[i] = N - I[i]

# Analytical solution, using exponential equation
Ianal = (K*I0*np.exp(r*t)) / (K+I0*(np.exp(r*t)-1))

plt.figure(0)
plt.plot(t, S, label='S')
plt.plot(t, I, label='I')
plt.plot(t, Ianal, label='I ana')
plt.legend(loc='best')
plt.xlabel('Time (hours)')
plt.ylabel('Population')
plt.savefig('q5-6_infectives_d=%d.png' % d, dpi=300, bbox_inches='tight')
plt.show()