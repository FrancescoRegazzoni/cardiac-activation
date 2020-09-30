#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

import matplotlib.pyplot as plt

def postprocess(output, filename):
    fig, axes = plt.subplots(2, 2, figsize=(6, 4))

    axes[0,0].plot(output['t'], output['Ca'])
    axes[0,0].set_ylabel('$[\mathrm{Ca}^{2+}]_i \, [\mu M]$')
    axes[0,0].set_ylim((0, 1.2))

    axes[0,1].plot(output['t'], output['SL'])
    axes[0,1].set_ylabel('$SL \, [\mu m]$')
    axes[0,1].set_ylim((2.0, 2.3))

    axes[1,0].plot(output['t'], output['P'])
    axes[1,0].set_ylabel('permissivity [-]')
    axes[1,0].set_ylim((0, 1))

    axes[1,1].plot(output['t'], output['Ta'])
    axes[1,1].set_ylabel('active tension [kPa]')
    axes[1,1].set_ylim((0, 80))

    for i in range(2):
        for j in range(2):
            axes[i,j].set_xlabel('t [s]')
            axes[i,j].set_xlim((output['t'][0], output['t'][-1]))

    fig.tight_layout()
    fig.savefig(filename)