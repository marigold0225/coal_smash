import matplotlib.pyplot as plt
import numpy as np


# Set global Matplotlib parameters
def set_matplotlib_params():
    plt.rcParams.update({
        'font.family': 'Times New Roman',
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.labelsize': 20,
        'ytick.labelsize': 20
    })


# Define particle configuration
particle_configs = {
    'p': {
        'xlim': [0.2, 2.5],
        'ylim': [10 ** (-7), 100],
        'file_prefix': 'p_pt_',
        'file_exp': '../tem/p.dat',
        'label': r'Proton',
    },
    'd': {
        'xlim': [0, 4.2],
        'ylim': [10 ** (-12), 10],
        'file_prefix': 'd_pt_',
        'file_exp': '../tem/d.dat',
        'label': r'Deutron',
    },
    'he4': {
        'xlim': [0.4, 5.2],
        'ylim': [10 ** (-12), 1],
        'file_prefix': 'alpha3_pt_',
        'file_exp': '../tem/he4.dat',
        'label': r'$^4$He',
    },
    'be': {
        'xlim': [0.4, 5.2],
        'ylim': [10 ** (-13), 10 ** (-4)],
        'file_prefix': 'Be3_pt_',
        'file_exp': '../tem/be.dat',
        'label': r'$^8$Be',
    },
    # Add more particles as needed
}


# Read and plot particle data
def read_and_plot_particle(ax, centrality, particle_type, scale_factors, markers):
    config = particle_configs[particle_type]
    filepath = f'../data/50000/{centrality}/{config["file_prefix"]}{centrality}.dat'
    dy = 0.1  # Define dy if it's constant or pass it as a parameter

    with open(filepath, 'r') as file:
        lines = file.readlines()

    rapidity_ranges = []
    data = []

    for line in lines:
        if line.startswith("Rapidity range:"):
            rapidity_values = line.split(", ")[1].split(":")[1]
            print(rapidity_values)
            if rapidity_values != 0:
                rapidity_ranges.append(line.strip().split(": ")[1])
                data.append([])
        elif line.strip():
            pt, density = map(float, line.split())
            data[-1].append([pt, density])

    for i, rapidity_data in enumerate(data):
        rapidity_data = np.array(rapidity_data)
        x = rapidity_data[:, 0]
        y = rapidity_data[:, 1] / dy * scale_factors[i]
        ax.plot(x, y, marker=markers[i], linestyle='-', label=rapidity_ranges[i], markersize=8)

    ax.set_title(f'{centrality}%-{config["label"]}', fontsize=20)
    ax.set_xlabel('Pt (GeV/c)', fontsize=22)
    ax.set_yscale('log')
    ax.set_xlim(config['xlim'])
    ax.set_ylim(config['ylim'])
    ax.grid(True)


# Main plotting routine
def main():
    set_matplotlib_params()
    rapidity_range = ["-0.1<y<0.0", "-0.2<y<-0.1", "-0.3<y<-0.2", "-0.4<y<-0.3", "-0.5<y<-0.4"
        , "-0.6<y<-0.5", "-0.7<y<-0.6", "-0.8<y<-0.7", "-0.9<y<-0.8", "-1.0<y<-0.9"]
    scale_factors = [10 ** (-i) for i in range(len(rapidity_range))]
    markers = ['o', 's', 'D', '^', 'v'
        , '<', '>', 'p', '*', 'h']

    # Define the centrality ranges
    centrality_ranges = ["0-10", "10-20", "20-40", "40-80"]
    fig, axs = plt.subplots(1, len(centrality_ranges), figsize=(20, 6), sharey=True)
    fig.suptitle(r'Au+Au @ FXT $\sqrt{s_{NN}}=3$ GeV', fontsize=22)

    # Plot the data for each centrality
    for i, centrality in enumerate(centrality_ranges):
        read_and_plot_particle(axs[i], centrality, 'd', scale_factors, markers)  # Change 'd' to desired particle key

        if i == 0:
            axs[i].set_ylabel(r'$d^2N/(2\pi p_T dy dp_t)[(GeV/c)^{-2}]$', fontsize=22)

    # Adjust layout and show
    axs[-1].legend(loc='upper right')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()


if __name__ == '__main__':
    main()
