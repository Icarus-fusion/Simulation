import numpy as np
import matplotlib.pyplot as plt
import math
plt.close('all')

T = 1.5e8
k_b = 1.380649e-23
M_d = 3.343583719e-27
M_t = 5.007356744e-27
facteur_d = np.sqrt(2 / math.pi) * (M_d / (k_b * T))**(3/2)
facteur_t = np.sqrt(2 / math.pi) * (M_t / (k_b * T))**(3/2)

x = np.linspace(0,int(5e6),int(int(5e6) / 10000))
yd = facteur_d * x**2 * np.exp(-M_d * x**2 / (2 * k_b * T))
yt = facteur_t * x**2 * np.exp(-M_t * x**2 / (2 * k_b * T))

plt.plot(x, yd, color='red', label='Deutérium')
plt.plot(x, yt, color='green', label='Tritium')
plt.grid()
plt.title("Probabilité de vitesse à 150 millions de degré")
plt.xlabel('v [m/s]')
plt.legend()
plt.show()