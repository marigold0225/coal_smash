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
M_proton = 0.938
M_deutron = 1.875
M_alpha = 3.727

# experimental data
proton_yields = 43.02
deuteron_yields = 5.41
alpha_yields = 0.21448


def read_and_plot_yield(ax_in, filepath, markers_in, rapidity_ranges, M_particle, g_particle, particle_label,
                        add_legend, add_color):
    # check if the file exists
    if not os.path.isfile(filepath):
        print(f"File {filepath} does not exist!")
        return
    # Open the file and read the lines
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Process the file
    label_yield_pairs = []  # This will store pairs of (label, yield)
    for line in lines:
        if line.startswith("Rapidity range:"):
            parts = line.strip().split(", ")
            label = parts[0].split(": ")[1]  # Extracting the rapidity range label
            yield_value = float(parts[1].split(": ")[1])  # Extracting the yield
            label_yield_pairs.append((label, yield_value))

    # Now plot only the yields for the rapidity ranges we are interested in
    for label, yield_value in label_yield_pairs:
        if label in rapidity_ranges and yield_value != 0.0:  # Check if the label is one we want
            marker = markers_in[rapidity_ranges.index(label)]  # Get the corresponding marker

            print(f"yield_value: {yield_value}, rapidity: {label}")
            if add_legend:
                ax_in.plot(M_particle, yield_value / g_particle / dy, marker=marker, linestyle='', label=label,
                           markersize=8, color=add_color[rapidity_ranges.index(label)])
            else:
                ax_in.plot(M_particle, yield_value / g_particle / dy, marker=marker, linestyle='', markersize=8,
                           color=add_color[rapidity_ranges.index(label)])

    ax_in.text(M_particle, max(yld for _, yld in label_yield_pairs) / g_particle / dy, particle_label, fontsize=16,
               ha='center', va='bottom')


# Define markers and rapidity ranges for each point
colors = ['k', 'r', 'b', 'm', 'c']
markers = ['o', 's', 'D', '*', 'p']
rapidity_ranges = ["-0.1<y<0.0", "-0.3<y<-0.2", "-0.5<y<-0.4",
                   "-0.7<y<-0.6", "-0.9<y<-0.8", ]

# Select the first centrality range
centrality = "0-10"
# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle('Particle Yield for Different Rapidity Ranges in 0-10 Centrality')

# Read and plot proton data
fileProton = f'../data/50000/{centrality}/p_pt_{centrality}.dat'
read_and_plot_yield(ax, fileProton, markers, rapidity_ranges, M_proton, g_proton, 'Proton', True, colors)

# Read and plot deuteron data
fileDeutron = f'../data/50000/{centrality}/d_pt_{centrality}.dat'
read_and_plot_yield(ax, fileDeutron, markers, rapidity_ranges, M_deutron, g_deutron, 'Deuteron', False, colors)

fileAlpha = f'../data/50000/{centrality}/alpha_pt_{centrality}.dat'
read_and_plot_yield(ax, fileAlpha, markers, rapidity_ranges, M_alpha, g_alpha, 'Alpha', False, colors)

# Set plot properties
ax.set_title(f'Au+Au {centrality}% 3GeV', fontsize=20)
ax.set_xlabel('Mass (GeV/c²)', fontsize=22)
ax.set_xlim(0.5, 4)
ax.set_yscale('log')
ax.set_ylim(0.06, 200)
# add experimental data
ax.scatter(M_proton, proton_yields, color='g', marker='*', s=100)  # 星号标记，更大尺寸
ax.scatter(M_deutron, deuteron_yields, color='g', marker='*', s=100)
ax.scatter(M_proton * 4, alpha_yields, color='g', marker='*', s=100)

# 添加文本说明实验数据
ax.text(M_proton, proton_yields, 'Proton (Exp)', color='g', verticalalignment='bottom')
ax.text(M_deutron, deuteron_yields, 'Deuteron (Exp)', color='g', verticalalignment='bottom')
ax.text(M_proton * 4, alpha_yields, 'Alpha (Exp)', color='g', verticalalignment='bottom')
ax.set_ylabel('dN/dy/(2J+1)', fontsize=22)
ax.legend(title="Rapidity Ranges", loc='upper right')
plt.show()
