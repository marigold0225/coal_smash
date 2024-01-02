import os

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
dy = 0.1
g_proton = 2
g_deutron = 3
g_alpha = 1
g_be = 1
M_proton = 0.938
M_deutron = 1.875
M_alpha = 3.727
M_Be = 8.005

# experimental data
p1 = 41.9
p2 = 37.320069688613415
p3 = 34.070154321567
d1 = 5.372987464539156
d2 = 5.070327024327639
d3 = 4.628791575712421
a1 = 0.21603910732275478
a2 = 0.2055652213825122
a3 = 0.18458102373095392


def read_and_plot_yield(ax_in, filepath, markers_in, rapidity_range, M_particle, g_particle, particle_label,
                        add_legend, add_color):
    # check if the file exists
    if not os.path.isfile(filepath):
        print(f"File {filepath} does not exist!")
        return
    # Open the file and read the lines
    with open(filepath, 'r') as file:
        lines = file.readlines()

    label_yield_pairs = []
    for line in lines:
        if line.startswith("Rapidity range:"):
            parts = line.strip().split(", ")
            label = parts[0].split(": ")[1]
            yield_value = float(parts[1].split(": ")[1])
            label_yield_pairs.append((label, yield_value))

    for label, yield_value in label_yield_pairs:
        if label in rapidity_range and yield_value != 0.0:
            marker = markers_in[rapidity_range.index(label)]

            if add_legend:
                ax_in.plot(M_particle, yield_value / g_particle / dy, marker=marker, linestyle='', label=label,
                           markersize=8, color=add_color[rapidity_range.index(label)])
            else:
                ax_in.plot(M_particle, yield_value / g_particle / dy, marker=marker, linestyle='', markersize=8,
                           color=add_color[rapidity_range.index(label)])

    ax_in.text(M_particle, max(yld for _, yld in label_yield_pairs) / g_particle / dy * 2, particle_label,
               fontsize=16,
               ha='center', va='bottom')


colors = ['k', 'r', 'b', 'm', 'c']
markers = ['o', 's', 'D', '*', 'p']
rapidity_ranges = ["-0.1<y<0.0", "-0.3<y<-0.2", "-0.5<y<-0.4"]
# , "-0.7<y<-0.6", "-0.9<y<-0.8", ]

centrality = "0-10"
fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle('Particle Yield for Different Rapidity Ranges in 0-10 Centrality')

fileProton = f'../data/50000/{centrality}/p_pt_{centrality}.dat'
read_and_plot_yield(ax, fileProton, markers, rapidity_ranges, M_proton, g_proton, 'p', True, colors)

fileDeutron = f'../data/50000/{centrality}/d_pt_{centrality}.dat'
read_and_plot_yield(ax, fileDeutron, markers, rapidity_ranges, M_deutron, g_deutron, 'd', False, colors)

fileAlpha = f'../data/50000/{centrality}/alpha_pt_{centrality}.dat'
read_and_plot_yield(ax, fileAlpha, markers, rapidity_ranges, M_alpha, g_alpha, r'$^4$He', False, colors)

fileBe = f'../data/50000/{centrality}/Be_pt_{centrality}.dat'
read_and_plot_yield(ax, fileBe, markers, rapidity_ranges, M_Be, g_be, 'Be', False, colors)

# Set plot properties
ax.set_title(f'Au+Au {centrality}% 3GeV', fontsize=20)
ax.set_xlabel('Mass (GeV/c²)', fontsize=22)
ax.set_xlim(0.5, 9)
ax.set_yscale('log')
ax.set_ylim(1e-7, 200)
# add experimental data
ax.scatter(M_proton, p1, color='g', marker='*', s=100)
ax.scatter(M_proton, p2, color='g', marker='*', s=100)
ax.scatter(M_proton, p3, color='g', marker='*', s=100)
ax.scatter(M_deutron, d1, color='g', marker='*', s=100)
ax.scatter(M_deutron, d2, color='g', marker='*', s=100)
ax.scatter(M_deutron, d3, color='g', marker='*', s=100)
ax.scatter(M_proton * 4, a1, color='g', marker='*', s=100)
ax.scatter(M_proton * 4, a2, color='g', marker='*', s=100)
ax.scatter(M_proton * 4, a3, color='g', marker='*', s=100)
ax.text(M_proton, p1, 'Proton (Exp)', color='g', verticalalignment='bottom')
ax.text(M_deutron, d1, 'Deuteron (Exp)', color='g', verticalalignment='bottom')
ax.text(M_proton * 4, a1, 'Alpha (Exp)', color='g', verticalalignment='bottom')
ax.set_ylabel('dN/dy/(2J+1)', fontsize=22)
ax.legend(title="Rapidity Ranges", loc='upper right')
plt.show()
