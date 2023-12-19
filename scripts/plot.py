import matplotlib.pyplot as plt
import numpy as np

alpha_two = np.loadtxt('../data/1000/d_mix_spv0-10.dat')
alpha_four = np.loadtxt('../data/1000/d_mix_spv10-20.dat')
deutron = np.loadtxt('../data/1000/d_mix_spv20-40.dat')
d2 = np.loadtxt('../data/1000/d_mix_spv40-80.dat')
x1 = alpha_two[:, 0]
y1 = alpha_two[:, 1]
x2 = alpha_four[:, 0]
y2 = alpha_four[:, 1]
x3 = deutron[:, 0]
y3 = deutron[:, 1]
x4 = d2[:, 0]
y4 = d2[:, 1]
plt.rc('font', family='Times New Roman')
plt.figure(figsize=(10, 6))
plt.plot(x1, y1, marker='o', linestyle='-', label='he4-dd')
plt.plot(x2, y2, marker='o', linestyle='-', label='he4-ppnn')
plt.plot(x3, y3, marker='o', linestyle='-', label='deutron')
plt.plot(x4, y4, marker='o', linestyle='-', label='d2')
plt.title('Transverse momentum distribution')
plt.xlabel('Pt(Gev/c)')
plt.ylabel('probability density')
plt.xlim(0, 4)
plt.ylim(1e-18, 10)
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()
