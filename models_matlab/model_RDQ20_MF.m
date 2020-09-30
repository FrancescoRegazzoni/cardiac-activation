function output = model_RDQ20_MF(input, params)
%__________________________________________________________________________
%
% This code implements the mean-field ODE model (MF-ODE) for cardiomyocytes
% force generation presented in [1, 2].
%
% Inputs:
%  - input: struct containing the input data
%     - input.times: time instants [s]
%     - input.Ca:    intracellular calcium ions concentration [micro M]
%     - input.SL:    sarcomere length [micro m]
%  - params: struct containing the model parameters
%
% Outputs:
%  - output: struct containing the output data
%     - output.times: time instants [s]
%     - output.Ca:    intracellular calcium ions concentration [micro M]
%     - output.SL:    sarcomere length [micro m]
%     - output.P:     permissivity [-]
%     - output.Ta:    active tension [kPa]
%
%__________________________________________________________________________
%
% Author: Francesco Regazzoni - MOX, Politecnico di Milano
% Email:  francesco.regazzoni@polimi.it
% Date:   2020
%__________________________________________________________________________
%
% References:
%
% [1] F. Regazzoni "Mathematical modeling and Machine Learning for the
%     numerical simulation of cardiac electromechanics", PhD Thesis -
%     Politecnico di Milano (2020)
%     http://hdl.handle.net/10589/152617
% [2] F. Regazzoni, L. Dede', A. Quarteroni "Biophysically detailed
%     mathematical models of multiscale cardiac active mechanics",
%     PLOS Computational Biology (2020)
%     https://doi.org/10.1371/journal.pcbi.1008294
%
%__________________________________________________________________________

    %% Model parameters
    if nargin == 1
        params = '../params/params_RDQ20-MF_human_body-temperature.json';
    end
    if ischar(params)
        params = jsondecode(fileread(params));
    end
    LA       = params.geometry.LA;             % [micro m]
    LM       = params.geometry.LM;             % [micro m]
    LB       = params.geometry.LB;             % [micro m]
    SL0      = params.geometry.SL0;            % [micro m]
    Q        = params.RU_steady_state.Q;       % [-]
    Kd0      = params.RU_steady_state.Kd0;     % [micro M]
    alphaKd  = params.RU_steady_state.alphaKd; % [micro M / micro m]
    mu       = params.RU_steady_state.mu;      % [-]
    gamma    = params.RU_steady_state.gamma;   % [-]
    Koff     = params.RU_kinetics.Koff;        % [s^-1]
    Kbasic   = params.RU_kinetics.Kbasic;      % [s^-1]
    r0       = params.XB_cycling.r0;           % [s^-1]
    alpha    = params.XB_cycling.alpha;        % [-]
    mu0_fP   = params.XB_cycling.mu0_fP;       % [s^-1]
    mu1_fP   = params.XB_cycling.mu1_fP;       % [s^-1]
    a_XB     = params.upscaling.a_XB;          % [kPa]

    %% Numerical parameters
    dt_RU = 2.5e-5;                                % [s]
    freqXB = round(1e-3 / dt_RU);                  % [-]
    freqenceRatesUpdate = round(2.5e-4 / 2.5e-5);  % [-]

    %% Initialization
    LMh = (LM - LB) / 2;
    frac_SO =  @(SL)((SL>LA)     .*(SL<=LM)     .*(SL-LA) + ...
                     (SL>LM)     .*(SL<=2*LA-LB).*(SL+LM-2*LA)*.5 + ...
                     (SL>2*LA-LB).*(SL<=2*LA+LB).*LMh + ...
                     (SL>2*LA+LB).*(SL<=2*LA+LM).*(LM+2*LA-SL)*.5 ...
               )/LMh;
    Kd = @(SL) Kd0 - alphaKd*(2.15-SL);
    kC = zeros(2,2);     %kC(A,a) = P( C_i^t+dt =~ A | (C_i,T_i)^t = (A,a) ) / dt
    kT = zeros(2,2,2,2); %kT(a,b,c,B) = P( T_i^t+dt =~ b | (T_i-1,T_i,T_i+1,C_i)^t = (a,b,c,B) ) / dt
    expMat = permute([0 1; 1 2],[1,3,2]);
    kC(2,1) = Koff;
    kC(2,2) = Koff/mu;
    kT(:,2,:,1) = Kbasic*gamma.^(2-expMat);
    kT(:,2,:,2) = Kbasic*gamma.^(2-expMat);
    kT(:,1,:,1) = Q*Kbasic/mu*gamma.^expMat;
    kT(:,1,:,2) = Q*Kbasic*gamma.^expMat;

    %% Time loop
    times = min(input.times):dt_RU:max(input.times);
    Ca    = interp1(input.times,input.Ca,times);
    SL    = interp1(input.times,input.SL,times);
    nT    = length(times);
    P     = zeros(1,nT);
    Ta    = zeros(1,nT);

    time_init = tic();
    fprintf('RDQ20-MF model. Computing... ')
    for iT = 1:nT
        if mod(iT-1,freqenceRatesUpdate) == 0
            update_rates_RU(Ca(iT),SL(iT));
        end
        if iT == 1
            x_RU = zeros(2,2,2,2); %x_RU(a,b,c,B) = P( (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) )
            x_RU(1,1,1,1) = 1;
            x_XB = zeros(2,2);
        else
            x_RU = x_RU + dt_RU*RU_rhs(x_RU);
            if mod(iT-1,freqXB) == 0
                dSLdt = (SL(iT)-SL(iT-freqXB))/(dt_RU*freqXB);
                x_XB = XB_advance(x_RU,x_XB,dt_RU*freqXB,dSLdt);
            end
        end
        P(iT) = sum(reshape(x_RU(:,2,:),[],1))*frac_SO(SL(iT));
        Ta(iT) = a_XB*sum(x_XB(2,:))*frac_SO(SL(iT));
    end
    fprintf('done. Time elapsed: %1.3f s\n', toc(time_init))

    output.t  = times';
    output.Ca = Ca';
    output.SL = SL';
    output.P  = P';
    output.Ta = Ta';

    %% Functions definition
    function update_rates_RU(Ca,SL)
        Kon = Koff/Kd(SL); % [micro M^-1 * s^-1]
        kC(1,:) = Kon*Ca;
    end

    function rhs = RU_rhs(x_RU)
        PhiT_C = x_RU.*kT;                    % ...(a,b,c,B) =  P( T_i^t+dt =~ b     , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        kT_L = permute(sum(sum(PhiT_C,1),4)./sum(sum(x_RU,1),4),[2,3,1,4]);  % ...(a,b,_,_) ~= P( T_{i-1}^t+dt =~ a | (T_{i-1},T_i)^t = (a,b) ) / dt
        kT_R = permute(sum(sum(PhiT_C,3),4)./sum(sum(x_RU,3),4),[3,1,2,4]);  % ...(_,b,c,_) ~= P( T_{i+1}^t+dt =~ c | (T_i,T_{i+1})^t = (b,c) ) / dt
        PhiT_L = kT_L.*x_RU;                  % ...(a,b,c,B) ~= P( T_{i-1}^t+dt =~ a , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        PhiT_R = kT_R.*x_RU;                  % ...(a,b,c,B) ~= P( T_{i+1}^t+dt =~ c , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        PhiC_C = x_RU.*permute(kC,[3,2,4,1]); % ...(a,b,c,B) =  P( C_i^t+dt =~ B     , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        PhiT_L(isnan(PhiT_L))=0;
        PhiT_L(isinf(PhiT_L))=0;
        PhiT_R(isnan(PhiT_R))=0;
        PhiT_R(isinf(PhiT_R))=0;
        rhs = - PhiT_L + flip(PhiT_L,1)...
              - PhiT_C + flip(PhiT_C,2)...
              - PhiT_R + flip(PhiT_R,3)...
              - PhiC_C + flip(PhiC_C,4);
    end

    function x_RU_new = XB_advance(x_RU, x_XB, dt, dSLdt)
        v = -dSLdt/SL0; % [1/s]
        perm = sum(reshape(x_RU(:,2,:),[],1));
        k_PN = sum(reshape(kT(:,2,:,:).*x_RU(:,2,:,:),[],1))./sum(reshape(x_RU(:,2,:,:),[],1));
        k_NP = sum(reshape(kT(:,1,:,:).*x_RU(:,1,:,:),[],1))./sum(reshape(x_RU(:,1,:,:),[],1));
        k_PN(isnan(k_PN)) = 0; k_NP(isnan(k_NP)) = 0; % if denominator is zero, also numerator is zero
        r = r0 + alpha*abs(v);
        diag_P = r + k_PN;
        diag_N = r + k_NP;
        sol = reshape(x_XB(:,:),4,1);
        rhs = [perm*mu0_fP; perm*mu1_fP; 0; 0];
        A = [-diag_P, 0, k_NP, 0; ...
             -v, -diag_P, 0, k_NP;...
             k_PN, 0, -diag_N, 0; ...
             0, k_PN, -v, -diag_N];
        sol_inf = -A\rhs;
        x_RU_new = reshape(sol_inf + expm(dt*A)*(sol - sol_inf),[2,2]);
    end

end