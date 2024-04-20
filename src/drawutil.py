import matplotlib.pyplot as plt


# https://brand.utexas.edu/identity/color/
utcolors = [
    '#f8971f', # orange
    '#ffd600', # yellow
    '#a6cd57', # light green
    '#579d42', # vivid green
    '#00a9b7', # teal
    '#005f86', # dark blue
    '#9cadb7', # blue gray
    '#d6d2c4', # tan gray
]

# resilience/flood level colors
rfcolors = [
    [0.659, 0.902, 0.114],
    [1.000, 0.949, 0.000],
    [1.000, 0.494, 0.000],
    [0.929, 0.110, 0.141],
    [0.328, 0.074, 0.531],
]


def mpl_config(scale):
    plt.rcParams['font.size'] = 9 * scale
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.labelsize'] = 9 * scale
    plt.rcParams['axes.titlesize'] = 9 * scale
    plt.rcParams['xtick.labelsize'] = 8 * scale
    plt.rcParams['ytick.labelsize'] = 8 * scale
    plt.rcParams['legend.title_fontsize'] = 9 * scale
    plt.rcParams['legend.fontsize'] = 8 * scale
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    # plt.rcParams['font.family'] = 'Times New Roman'
    # mpl.use("pgf")
    # pgf_with_pdflatex = {
    #     "pgf.texsystem": "pdflatex",
    #     "pgf.preamble": r"""
    #     \usepackage{amsmath}
    #     """
    # }
    # mpl.rcParams.update(pgf_with_pdflatex)


def display_utcolors():
    plt.bar(range(len(utcolors)), [1] * len(utcolors), color=utcolors)
    plt.show()


def display_rfcolors():
    plt.bar(range(len(rfcolors)), [1] * len(rfcolors), color=rfcolors)
    plt.show()
