% Author: Francesco Regazzoni - MOX, Politecnico di Milano
% Email:  francesco.regazzoni@polimi.it
% Date:   2020

clear

%% Time interval
Tmax = .6;   % [s]

%% Calcium transient
c0 = .1;     % [micro M]
cmax = .9;  % [micro M]
tau1 = .02;  % [s]
tau2 = .05;  % [s]
t0 = 0.01;   % [s]
beta = (tau1/tau2)^(-1/(tau1/tau2 - 1)) - (tau1/tau2)^(-1/(1 - tau2/tau1));
Ca_base = @(t) c0 + (t>=t0) .* ((cmax - c0) / beta * (exp(-(t-t0)/tau1) - exp(-(t-t0)/tau2)));

%% SL transient
SL0 = 2.2;     % [micro m]
SL1 = SL0*.97; % [micro m]
SLt0 = .05 ;   % [s]
SLt1 = .35;    % [s]
SLtau0 = .05;  % [s]
SLtau1 = .02;  % [s]
SL_base = @(t) SL0+ (SL1-SL0) * (max(0,1-exp((SLt0-t)/SLtau0)) - max(0,1-exp((SLt1-t)/SLtau1)));

%% Input definition
input.times = 0:1e-4:Tmax;
input.Ca = Ca_base(input.times);
input.SL = SL_base(input.times);

%% Simulation (select here the model)
% output = model_RDQ18(input);
% output = model_RDQ20_SE(input);
output = model_RDQ20_MF(input);

%% Postprocessing
postprocess(output)