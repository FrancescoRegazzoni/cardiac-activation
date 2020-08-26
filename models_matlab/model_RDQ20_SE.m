function output = model_RDQ20_SE(input, params)
%__________________________________________________________________________
%
% This code implements the spatially-explicitly ODE model (SE-ODE) for
% cardiomyocytes force generation presented in [1, 2].
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
%     submitted (2020)
%     https://arxiv.org/abs/2004.07910
%
%__________________________________________________________________________

    %% Model parameters
    if nargin == 1
        params = '../params/params_RDQ20-SE_human_body-temperature.json';
    end
    if ischar(params)
        params = jsondecode(fileread(params));
    end
    n_RU     = params.geometry.n_RU;           % [-]
    LA       = params.geometry.LA;             % [micro m]
    LM       = params.geometry.LM;             % [micro m]
    LB       = params.geometry.LB;             % [micro m]
    SL0      = params.geometry.SL0;            % [micro m]
    Lsmooth  = params.geometry.Lsmooth;        % [micro m]
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
    y_j = @(j) LA * (j-.5)/n_RU;
    yLA = @(SL) 2*LA - SL;
    yM0 = @(SL) (2*LA - SL + LB) / 2;
    yM1 = @(SL) (2*LA - SL + LM) / 2;
    ChiMF=@(SL,i) .5*tanh((y_j(i) - yM0(SL))/Lsmooth) + .5*tanh(-(y_j(i) - yM1(SL))/Lsmooth);
    ChiSF=@(SL,i) .5 + .5*tanh((y_j(i) - yLA(SL))/Lsmooth);
    Kd = @(SL) Kd0 - alphaKd*(2.15-SL);

    Kpn0 = Kbasic*gamma*gamma;
    Kpn1 = Kbasic*gamma*gamma;
    Knp0 = Q*Kbasic/mu;
    Knp1 = Q*Kbasic;
    K1off= Koff/mu;
    ratesRU       = zeros(n_RU,4,4,4,4); %ratesRU(i,a,b,c,D) = P((a,D,c)_i^t+dt | (a,b,c)_i^t ) / dt
    ones_414      = ones(n_RU,4,1,4);
    expMat        = repmat(permute([0 0 1 1; 0 0 1 1; 1 1 2 2; 1 1 2 2],[3,1,4,2]),n_RU,1,1,1);
    jMat          = repmat((1:n_RU)',1,4,1,4);
    gammaPlusExp  = gamma.^expMat;
    gammaMinusExp = gamma.^-expMat;
    ratesRU(:,:,2,:,1) = Koff*ones_414;
    ratesRU(:,:,3,:,4) = K1off*ones_414;
    ratesRU(:,:,4,:,1) = Kpn0*gammaMinusExp;
    ratesRU(:,:,3,:,2) = Kpn1*gammaMinusExp;
    P_local = @(x_RU)  [sum(sum(sum(x_RU(1,3:4,:,:),2),3),4); ... % first unit
                        sum(sum(sum(x_RU(:,:,3:4,:),2),3),4); ... % central units
                        sum(sum(sum(x_RU(end,:,:,3:4),2),3),4)];  % last unit

    %% Time loop
    times = min(input.times):dt_RU:max(input.times);
    Ca    = interp1(input.times,input.Ca,times);
    SL    = interp1(input.times,input.SL,times);
    nT    = length(times);
    P     = zeros(1,nT);
    Ta    = zeros(1,nT);

    time_init = tic();
    fprintf('RDQ20-SE model. Computing... ')
    for iT = 1:nT
        if mod(iT-1,freqenceRatesUpdate) == 0
            update_rates_RU(Ca(iT),SL(iT));
        end
        if iT == 1
            x_RU = zeros(n_RU-2,4,4,4); %x_RU(i,a,b,c) = P((X_{i-1},X_i,X_{i+1}) = (a,b,c) )
            x_RU(:,1,1,1) = 1;
            x_XB = zeros(n_RU,2,2);
        else
            x_RU = x_RU + dt_RU*RU_rhs(x_RU);
            if mod(iT-1,freqXB) == 0
                dSLdt = (SL(iT)-SL(iT-freqXB))/(dt_RU*freqXB);
                x_XB = XB_advance(x_RU,x_XB,dt_RU*freqXB,SL(iT),dSLdt);
            end
        end
        P(iT) = mean(P_local(x_RU));
        Ta(iT) = a_XB*mean(sum(x_XB(:,2,:),3));
    end
    fprintf('done. Time elapsed: %1.3f s\n', toc(time_init))

    output.times = times;
    output.Ca    = Ca;
    output.SL    = SL;
    output.P     = P;
    output.Ta    = Ta;

    %% Functions definition
    function update_rates_RU(Ca,SL)
        Kon = Koff/Kd(SL); % [micro M^-1 * s^-1]
        Chi_2_gammaPlus = ChiSF(SL,jMat).*gammaPlusExp;
        ratesRU(:,:,1,:,2) = Kon*Ca;
        ratesRU(:,:,4,:,3) = Kon*Ca;
        ratesRU(:,:,2,:,3) = Knp1*Chi_2_gammaPlus;
        ratesRU(:,:,1,:,4) = Knp0*Chi_2_gammaPlus;
    end

    function rhs = RU_rhs(x_RU)
        x2_L = sum(x_RU(2:end,:,:,:),4);                       % xODE2_L(i-2,a,b) = P( (X_{i-1},X_i) = (a,b) ) , i = 3:RU_nu-1
        x2_R = permute(sum(x_RU(1:end-1,:,:,:),2),[1 3 4 2]);  % xODE2_R(i-1,b,c) = P( (X_i,X_{i+1}) = (b,c) ) , i = 2:RU_nu-2
        PhiC = ratesRU(2:(n_RU-1),:,:,:,:) .* x_RU;
        PsiL_squeezed(2:(n_RU-2),:,:,:,:) = permute(sum(PhiC(1:(n_RU-3),:,:,:,:),2),[1,3,4,6,5,2])./ x2_L;
        PsiL_squeezed(isnan(PsiL_squeezed))=0;
        PsiL_squeezed(isinf(PsiL_squeezed))=0;
        PsiL_squeezed(1,:,:,:,:) = permute(ratesRU(1,1,:,:,:),[1,3,4,6,5,2]);
        PhiL = PsiL_squeezed .* x_RU;
        PsiR_squeezed(1:(n_RU-3),:,:,:,:) = permute(sum(PhiC(2:(n_RU-2),:,:,:,:),4),[1,6,2,3,5,4]) ./ permute(x2_R,[1,4,2,3]);
        PsiR_squeezed(isnan(PsiR_squeezed))=0;
        PsiR_squeezed(isinf(PsiR_squeezed))=0;
        PsiR_squeezed(n_RU-2,:,:,:,:) = permute(ratesRU(n_RU,:,:,1,:),[1,6,2,3,5,4]);
        PhiR = PsiR_squeezed .* x_RU;
        rhs = sum( permute(PhiC,[1,2,5,4,3]) - PhiC + permute(PhiL,[1,5,3,4,2]) - PhiL + permute(PhiR,[1,2,3,5,4]) - PhiR,5);
    end

    function x_XB_new = XB_advance(x_RU, x_XB, dt, SL, dSLdt)
        v = -dSLdt/SL0; % [1/s]
        x_XB_new = zeros(n_RU,2,2);

        perm = P_local(x_RU);
        for i = 1:2
            if i == 1
                a1 = 3; a2 = 4; b1 = 2; b2 = 1;
                a_marginal = perm;
            else
                a1 = 2; a2 = 1; b1 = 3; b2 = 4;
                a_marginal = 1 - perm;
            end
            k_perm = [ sum(sum(permute(ratesRU(1,1,a1,:,b1),[4,1,2,3,5]) .*  permute(x_RU(1,a1,:,:),[3,4,1,2]),1),2) + ...
                       sum(sum(permute(ratesRU(1,1,a2,:,b2),[4,1,2,3,5]) .*  permute(x_RU(1,a2,:,:),[3,4,1,2]),1),2) ; ...
                       sum(sum(ratesRU(2:end-1,:,a1,:,b1).*x_RU(:,:,a1,:),2),4) + ...
                       sum(sum(ratesRU(2:end-1,:,a2,:,b2).*x_RU(:,:,a2,:),2),4) ; ...
                       sum(sum(permute(ratesRU(end,:,a1,1,b1),[2,1,3,4,5]) .*  permute(x_RU(end,:,:,a1),[3,2,1,4]),1),2) + ...
                       sum(sum(permute(ratesRU(end,:,a2,1,b2),[2,1,3,4,5]) .*  permute(x_RU(end,:,:,a2),[3,2,1,4]),1),2) ];
            k_perm = k_perm ./ a_marginal;
            if i == 1
                k_PN = k_perm;
            else
                k_NP = k_perm;
            end
        end
        chiSF_i = ChiSF(SL,(1:n_RU)').*ChiMF(SL,(1:n_RU)');
        mu0_fP_i = mu0_fP*chiSF_i;
        mu1_fP_i = mu1_fP*chiSF_i;
        r = r0 + alpha*abs(v);
        diag_P = r + k_PN;
        diag_N = r + k_NP;

        for i = 1:n_RU
            sol_i = reshape(x_XB(i,:,:),4,1);
            rhs = [perm(i)*mu0_fP_i(i); perm(i)*mu1_fP_i(i); 0; 0];
            A = [-diag_P(i), 0, k_NP(i), 0; ...
                 -v, -diag_P(i), 0, k_NP(i);...
                 k_PN(i), 0, -diag_N(i), 0; ...
                 0, k_PN(i), -v, -diag_N(i)];
            sol_inf_i = -A\rhs;
            x_XB_new(i,:,:) = reshape(sol_inf_i + expm(dt*A)*(sol_i - sol_inf_i),[2,2]);
        end
    end

end