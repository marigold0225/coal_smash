import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def set_matplotlib_params():
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.labelsize'] = 20
    matplotlib.rcParams['ytick.labelsize'] = 20


def read_and_plot_pt(ax, filepath, scale_factor, mark, color_in, exp, dy):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    rapidity_ranges = []
    data = []

    for line in lines:
        if line.startswith("Rapidity range:"):
            range_parts = line.split(",")[0].split(":")[1]
            rapidity_ranges.append(range_parts)
            data.append([])
        elif line.strip():
            pt, density = map(float, line.split())
            data[-1].append([pt, density])

    for i, (rapidity_data, scale, marker, c) in enumerate(zip(data, scale_factor, mark, color_in)):
        rapidity_data = np.array(rapidity_data)
        if exp:
            x = rapidity_data[:, 0]
            y = rapidity_data[:, 1] / dy * scale
            ax.plot(x, y, marker=marker, linestyle='--', markersize=8, color=c)
        else:
            x = rapidity_data[:, 0]
            y = rapidity_data[:, 1] / dy * scale
            ax.plot(x, y, marker=marker, linestyle='-', label=rapidity_ranges[i], markersize=8, color=c)


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
        'ylim': [10 ** (-9), 10],
        'file_prefix': 'd_pt_',
        'file_exp': '../tem/d.dat',
        'label': r'Deutron',
    },
    'he4': {
        'xlim': [0.4, 5.2],
        'ylim': [10 ** (-10), 1],
        'file_prefix': 'alpha_pt_',
        'file_exp': '../tem/he4.dat',
        'label': r'$^4$He',
    },
    'be': {
        'xlim': [0.4, 5.2],
        'ylim': [10 ** (-13), 10 ** (-4)],
        'file_prefix': 'Be_pt_',
        'file_exp': '../tem/be.dat',
        'label': r'$^8$Be',
    },
}


def plot_particle(ax, particle_type, centrality, scale_factors, markers, colors):
    config = particle_configs[particle_type]
    filedata = f'../data/50000/{centrality}/{config["file_prefix"]}{centrality}.dat'
    read_and_plot_pt(ax, filedata, scale_factors, markers, colors, False, dy=0.1)
    filedata_exp = f'{config["file_exp"]}'
    scale_factor_exp = [1, 1, 1, 1, 1]
    markers_exp = ['x', 'x', 'x', 'x', 'x']
    read_and_plot_pt(ax, filedata_exp, scale_factor_exp, markers_exp, colors, True, dy=1)
    ax.set_title(f'{centrality}%-{config["label"]}', fontsize=20)
    ax.set_xlabel('Pt (GeV/c)', fontsize=22)
    ax.set_yscale('log')
    ax.set_xlim(config['xlim'])
    ax.set_ylim(config['ylim'])
    ax.grid(True)
    ax.set_ylabel(r'$d^2N/(2\pi p_T dy dp_t)[(GeV/c)^{-2}]$', fontsize=22)
    ax.legend(loc='upper right')


def main():
    set_matplotlib_params()
    fig, ax = plt.subplots(figsize=(10, 6))  # particle type: p, d, he4, be
    plot_particle(ax, 'he4', '0-10', scale_factors=[10 ** (-i) for i in range(5)],
                  markers=['o', 's', 'D', '^', 'v'], colors=['k', 'r', 'b', 'm', 'c'])
    plt.show()


if __name__ == '__main__':
    main()
