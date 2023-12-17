import matplotlib.pyplot as plt
import numpy as np

alpha_two = np.loadtxt('../tem/alpha_mix_spv_two.dat')
# alpha_four = np.loadtxt('../tem/100/alpha_mix_spv_four.dat')
deutron = np.loadtxt('../tem/d_mix_spv.dat')
x1 = alpha_two[:, 0]
y1 = alpha_two[:, 1]
# x2 = alpha_four[:, 0]
# y2 = alpha_four[:, 1]
x3 = deutron[:, 0]
y3 = deutron[:, 1]
plt.rc('font', family='Times New Roman')
plt.figure(figsize=(10, 6))
plt.plot(x1, y1, marker='o', linestyle='-', label='two-component')
# plt.plot(x2, y2, marker='o', linestyle='-', label='four-component')
plt.plot(x3, y3, marker='o', linestyle='-', label='deutron')
plt.title('Transverse momentum distribution')
plt.xlabel('Pt(Gev/c)')
plt.ylabel('probability density')
plt.xlim(0, 4)
plt.ylim(1e-18, 10)
plt.yscale('log')
plt.grid(True)
plt.show()
