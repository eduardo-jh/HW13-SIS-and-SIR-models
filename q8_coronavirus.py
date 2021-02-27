#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW13 - Question 8. SIR model for Coronavirus

Created on Sat Feb 27 14:40:30 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

d = 24*14  # time of infection, hours
alpha = 1.5e-4  #rate of infection 1/person-hr
N = 500  # total population
B = 1/d  # rate of recovery
steps = 550  # time of simulation
dt = 1  # time step, hours
I0 = 1
S0 = N - I0

t = np.linspace(0, steps, int(steps/dt)+1)

S = np.zeros(len(t))
I = np.zeros(len(t))
R = np.zeros(len(t))
S[0], I[0], R[0] = S0, I0, 0

# Numerical solution to SIR model, using Euler method
for i in range(1, len(t)):
    S[i] = S[i-1] + (-alpha * S[i-1] * I[i-1]) * dt
    I[i] = I[i-1] + (alpha * S[i-1] * I[i-1] - B*I[i-1]) * dt
    R[i] = R[i-1] + B * I[i-1] * dt

plt.figure(0)
plt.plot(t, S, 'b-', label='S')
plt.plot(t, I, 'r--', label='I')
plt.plot(t, R, 'g:', label='R')
plt.legend(loc='best')
plt.xlabel('Time (hours)')
plt.ylabel('Population')
plt.savefig('q8_SIR_coronavirus.png', dpi=300, bbox_inches='tight')
plt.show()