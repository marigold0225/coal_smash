import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('../tem/alpha_mix_spv.dat')

x = data[:, 0]
y = data[:, 1]
plt.rc('font', family='Times New Roman')
plt.figure(figsize=(10, 6))
plt.plot(x, y, 'r-', label='d_mix_spv')
plt.title('Transverse momentum distribution')
plt.xlabel('Pt(Gev/c)')
plt.ylabel('probability density')
plt.xlim(0, 4)
plt.yscale('log')
plt.grid(True)
plt.show()
