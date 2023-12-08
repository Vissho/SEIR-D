#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def draw(filenames, labels):

    plt.rcParams["legend.markerscale"] = 2.0
    plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = ['PT Sans']
    plt.rcParams['font.size'] = '12'
    plt.rcParams["legend.loc"] = "upper left" # 'upper left', 'upper right', 'lower left', 'lower right', 'center'
    plt.rcParams["legend.fontsize"] = "8"
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(10*cm, 7*cm))
    ax = fig.add_subplot(111)
    ax.set_title("")
    ax.set(xlabel="Выявленные случаи (чел)", ylabel="Время (д)")
    ax.label_outer()
    # ax.set_xlim(right=110)
    # ax.xaxis.set_ticks(np.arange(min(x), max(x)+1, 1))
    # ax.yaxis.set_ticks(np.arange(min(y), max(y)+1, 1))
    ax.set_ylim([0, 104])
    ax.set_yticks(np.arange(0, ax.get_ylim()[1]+1, 15))
    ax.xaxis.set_tick_params(direction='in', which='both')
    ax.yaxis.set_tick_params(direction='in', which='both')
    #ax.loglog()

    i = 0
    for (fname, datalabel) in zip(filenames, labels):
        data = np.loadtxt(fname, ndmin=2)
        x = data[:, 1]
        y = data[:, 0]

        if datalabel == "xi":
            ax.scatter(x, y, c="red", s=5,label=datalabel)
        else:
            ax.plot(x,y,'-o',markersize=1,linewidth=0.5,label=datalabel)

    plt.tight_layout()
    ax.legend()
    fig.savefig('seir-d.png', dpi=300)
    #fig.savefig('chart.pdf', dpi=300)

if __name__ == "__main__":
    draw(["data.dat"], ["23 марта - 31 мая 2020 г."])