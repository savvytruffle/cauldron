import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

def user_rc(lw=1.5):
    """Set plotting RC parameters"""
    # These are the "Tableau 20" colors as RGB.
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    # Change some of the default line-widths, font sizes for xticks, labels, titles,
    # and the color cycle to tableau20
    plt.rc('lines', linewidth=lw)
    plt.rc('font', size=14, weight='normal')
    plt.rc('xtick', labelsize=14)
    plt.rc('xtick.major', size=6, width=1)
    plt.rc('axes', prop_cycle=cycler(c=tableau20), lw=1, labelsize=18, titlesize=22)
    plt.rc('figure', titlesize=22, figsize=(10,8))
    return tableau20


# storing tableau20 color cycle into colors array, in case you want to use specific ones
# i.e., plt.plot(np.arange(10), color=colors[5])
colors=user_rc()

def example_plot():
    x = np.linspace(0, 2*np.pi, 100)
    y1 = np.sin(x)
    y2 = np.sin(x+0.5)
    y3 = 2.0*np.sin(x+0.5)
    y4 = 2.7*np.cos(x)
    
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)
    plt.plot(x, y4)
    plt.xlabel('x axis label here')
    plt.ylabel('y axis label here')
    plt.title('spaceman spiff strikes again!')
    plt.show()