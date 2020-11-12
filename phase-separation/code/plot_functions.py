import math
import numpy as np
import matplotlib.pyplot as plt


def plot_function(filename, graphname, title, x, y):

    X = np.loadtxt(filename)[:,0]
    Y = np.loadtxt(filename)[:,1]

    plt.plot(X, Y, linewidth = 0.5)
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()

def plot_function_Y1Y2_vs_X(y1_filename, y2_filename, graphname, title, x, y, y1, y2):
                            
    Y1 = np.loadtxt(y1_filename)[:,1]
    Y2 = np.loadtxt(y2_filename)[:,1]
    X = np.loadtxt(y1_filename)[:,0]

    plt.plot(X, Y1, linewidth = 0.5, label = y1)
    plt.plot(X, Y2, linewidth = 0.5, label = y2)
    plt.legend()
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()

def plot_function_Y1Y2_vs_X_different_scales(y1_filename, y2_filename, graphname, title, x, y1, y2):

    Y1 = np.loadtxt(y1_filename)[:,1]
    Y2 = np.loadtxt(y2_filename)[:,1]
    X = np.loadtxt(y1_filename)[:,0]

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel(x)
    ax1.set_ylabel(y1, color=color)
    ax1.set(ylim = (0, 80000))
    ax1.plot(X, Y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that 
                       #shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(y2, color=color)
    ax2.set(ylim = (0, 700))
    ax2.plot(X, Y2, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    plt.title(title)
    plt.savefig(graphname)
    plt.show()
        
def plot_function_w_errors(filename, graphname, title, x, y):

    X = np.loadtxt(filename)[:,0]
    Y = np.loadtxt(filename)[:,1]
    errors = np.loadtxt(filename)[:,2]

    plt.plot(X, Y, linewidth = 0.5)
    plt.errorbar(X, Y, linewidth = 0.5, yerr = errors, 
        ecolor = 'red')
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()

def plot_heatmap(filename, graphname, title, x, y):

    fig, ax = plt.subplots(1,1)
    data = np.loadtxt(filename)
    data = np.transpose(data)
    im=plt.imshow(data, origin = 'lower')
    plt.colorbar(im)
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()

def plot_quiver(x_filename, y_filename, graphname, title, x, y, n):

    xvalues = np.loadtxt(x_filename)
    yvalues = np.loadtxt(y_filename)

    ecks, why = np.meshgrid(np.arange(0, n, 1), np.arange(0, n ,1))

    plt.quiver(why, ecks, xvalues, yvalues, linewidth = 0.5)
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()