import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16


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
        y = rapidity_data[:, 1] / 0.2 * scale
        ax.plot(x, y, marker=marker, linestyle='-', label=rapidity_ranges[i], markersize=8)


# Define the scale factors and markers for each rapidity range
rapidity_range = ["-0.1<y<0.0", "-0.2<y<-0.1", "-0.3<y<-0.2", "-0.4<y<-0.3", "-0.5<y<-0.4", "-0.6<y<-0.5",
                  "-0.7<y<-0.6", "-0.8<y<-0.7", "-0.9<y<-0.8", "-1.0<y<-0.9"]
scale_factors = [10 ** (-i) for i in range(len(rapidity_range))]
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h']

# Define the centrality ranges
centrality_ranges = ["0-10", "10-20", "20-40", "40-80"]

# Create the subplot
fig, axs = plt.subplots(1, len(centrality_ranges), figsize=(20, 6), sharey=True)
fig.suptitle('Transverse momentum distribution for different rapidity ranges')

# Plot the data for each centrality
for i, centrality in enumerate(centrality_ranges):
    filename = f'../data/6/{centrality}/d_pt_{centrality}.dat'
    read_and_plot_pt(axs[i], filename, scale_factors, markers)

    axs[i].set_title(f'Centrality {centrality}')
    axs[i].set_xlabel('Pt (GeV/c)')
    axs[i].set_yscale('log')
    axs[i].set_xlim(0, 2.5)
    axs[i].set_ylim(10 ** (-15), 100)
    axs[i].grid(True)
    if i == 0:
        axs[i].set_ylabel('Probability Density (a.u.)')

# Adjust the layout and add a legend to the last subplot
# axs[-1].legend(loc='upper right')
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
