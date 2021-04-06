import numpy as np
import matplotlib.pyplot as plt

time, a, b, c = np.loadtxt('data.txt', skiprows=1).T

plt.plot(time, a, label='species A', c='green')
plt.plot(time, b, label='species B', c='red')
plt.plot(time, c, label='C', c='blue')
#plt.plot(time, d, label='D', c='green')
#plt.plot(time, e, label='E', c='yellow')
plt.legend(loc='upper left')
plt.xlabel('Time')
plt.ylabel('Amounts')
plt.title('Gillespie Algorithm - Stochastic Simulation of Chemical Reactions')
plt.grid()
plt.show()
