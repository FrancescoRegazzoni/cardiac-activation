#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

from model_RDQ18 import model_RDQ18
from model_RDQ20_SE import model_RDQ20_SE
from model_RDQ20_MF import model_RDQ20_MF
from postprocess import postprocess
import numpy as np

#%% Time interval
Tmax = .6   # [s]

#%% Calcium transient
c0 = .1     # [micro M]
cmax = 0.9  # [micro M]
tau1 = .02  # [s]
tau2 = .05  # [s]
t0 = 0.01   # [s]
beta = (tau1/tau2)**(-1/(tau1/tau2 - 1)) - (tau1/tau2)**(-1/(1 - tau2/tau1))
Ca_base = lambda t: c0 + (t>=t0) * ((cmax - c0) / beta * (np.exp(-(t-t0)/tau1) - np.exp(-(t-t0)/tau2)))

#%% SL transient
SL0 = 2.2;     # [micro m]
SL1 = SL0*.97; # [micro m]
SLt0 = .05;    # [s]
SLt1 = .35;    # [s]
SLtau0 = .05;  # [s]
SLtau1 = .02;  # [s]
SL_base = lambda t: SL0 + (SL1-SL0) * (np.clip(1-np.exp((SLt0-t)/SLtau0),0,None) - np.clip(1-np.exp((SLt1-t)/SLtau1),0,None))

#%% Input definition
inputs = {}
inputs['times'] = np.arange(0,Tmax,1e-4)
inputs['Ca'] = Ca_base(inputs['times']);
inputs['SL'] = SL_base(inputs['times']);

#%% Simulation (select here the model)
# model = model_RDQ18()
# model = model_RDQ20_SE()
model = model_RDQ20_MF()

output = model.solve(inputs)

#%% Postprocessing
postprocess(output)