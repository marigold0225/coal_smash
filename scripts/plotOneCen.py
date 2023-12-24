import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20

dy = 0.1
# experimental data path
p = '../tem/p.csv'
p1 = '../tem/p1.csv'
d = '../tem/deutron.csv'
d1 = '../tem/deutron1.csv'
he4 = '../tem/he4.csv'
he41 = '../tem/he41.csv'
# x and limit
p_xlim = [0.2, 2.5]
p_ylim = [10 ** (-15), 100]
d_xlim = [0, 4.2]
d_ylim = [10 ** (-7), 10]
he4_xlim = [0.4, 5.2]
he4_ylim = [10 ** (-8), 1]


# experimental data
def plot_exp_data(ax_in, filepath, marker='x', color='r'):
    data = pd.read_csv(filepath, header=None)
    data.columns = ['pt', 'density']
    x = data['pt'].values
    y = data['density'].values
    ax_in.plot(x, y, marker=marker, linestyle='', markersize=8, color=color, label='exp data')


def read_and_plot_pt(ax, filepath, scale_factors, markers):
    # Open the file and read the lines
    with open(filepath, 'r') as file:
        lines = file.readlines()

    rapidity_ranges = []
    data = []

    # Process the file
    for line in lines:
        if line.startswith("Rapidity range:"):
            # If a new rapidity range is found, reset the data list and store the range
            rapidity_ranges.append(line.strip().split(": ")[1])
            data.append([])
        elif line.strip():
            # If it's a data line, append the data to the last rapidity range's list
            pt, density = map(float, line.split())
            data[-1].append([pt, density])

    # Plot the data
    for i, (rapidity_data, scale, marker) in enumerate(zip(data, scale_factors, markers)):
        rapidity_data = np.array(rapidity_data)
        x = rapidity_data[:, 0]
        y = rapidity_data[:, 1] / dy * scale
        # ax.plot(x, y, marker=marker, linestyle='-', label=rapidity_ranges[i], markersize=8)
        ax.scatter(x, y, marker=marker, linestyle='-', label=rapidity_ranges[i], s=20)


# Define the scale factors and markers for each rapidity range
rapidity_range = ["-0.1<y<0.0", "-0.2<y<-0.1"]
scale_factors = [10 ** (-i) for i in range(len(rapidity_range))]
markers = ['o', 's', 'D', '^', 'v']

centrality_ranges = ["0-10"]

fig, ax = plt.subplots(figsize=(10, 6))  # 注意这里使用了ax而不是axs
fig.suptitle(r'Au+Au @ FXT $\sqrt{s_{NN}}=3$ GeV :he4', fontsize=22)

centrality = centrality_ranges[0]
filename = f'../data/50000/{centrality}/Be_pt_{centrality}.dat'  ##p_pt , alpha_pt, d_pt
read_and_plot_pt(ax, filename, scale_factors, markers)

ax.set_title(f'{centrality}%', fontsize=20)
ax.set_xlabel('Pt (GeV/c)', fontsize=22)
ax.set_yscale('log')
# ax.set_xlim(d_xlim)
# ax.set_ylim(d_ylim)
ax.grid(True)
ax.set_ylabel(r'$d^2N/(2\pi p_T dy dp_t)[(GeV/c)^{-2}]$', fontsize=22)

# 绘制实验数据
# plot_exp_data(ax, d, marker='x', color='r')
# plot_exp_data(ax, d1, marker='x', color='b')

# 添加图例
ax.legend(loc='upper right')
plt.show()
