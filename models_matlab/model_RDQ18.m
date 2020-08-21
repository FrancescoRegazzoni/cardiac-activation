function output = model_RDQ18(input, params)
%__________________________________________________________________________
%
% This code implements the ODE model for sarcomere activation proposed
% in [1].
%
% Inputs:
%  - input: struct containing the input data
%     - input.times: time instants [s]
%     - input.Ca:    intracellular calcium ions concentration [micro M]
%     - input.SL:    sarcomere length [micro m]
%  - params: struct containing the model parameters, or path of a json file
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
% [1] F. Regazzoni, L. Dedè, A. Quarteroni "Active contraction of cardiac
%     cells: a reduced model for sarcomere dynamics with cooperative
%     interactions", Biomechanics and Modeling in Mechanobiology (2018)
%     https://doi.org/10.1007/s10237-018-1049-0
%__________________________________________________________________________

    %% Model parameters
    if nargin == 1
        params = '../params/params_RDQ18.json';
    end
    if ischar(params)
        params = jsondecode(fileread(params));
    end
    LA       = params.geometry.LA;             % [micro m]
    LM       = params.geometry.LM;             % [micro m]
    LB       = params.geometry.LB;             % [micro m]
    n_RU     = params.geometry.n_RU;           % [-]
    Lsmooth  = params.geometry.Lsmooth;        % [micro m]
    Q0       = params.RU_steady_state.Q0;      % [-]
    SLQ      = params.RU_steady_state.SLQ;     % [micro m]
    alphaQ   = params.RU_steady_state.alphaQ;  % [micro m^-1]
    mu       = params.RU_steady_state.mu;      % [-]
    gamma    = params.RU_steady_state.gamma;   % [-]
    Kon      = params.RU_kinetics.Kon;         % [micro M^-1 * s^-1]
    Koff     = params.RU_kinetics.Koff;        % [s^-1]
    Kbasic   = params.RU_kinetics.Kbasic;      % [s^-1]
    TaMax    = params.upscaling.TaMax;         % [kPa]

    %% Numerical parameters
    dt = 2.5e-5;               % [s]
    freqenceRatesUpdate = 10;  % [-]

    %% Initialization
    x_i  = @(i) (LM-LB)*.5 * i/n_RU;
    xAZ  = @(SL) (SL-LB)/2;
    xLA  = @(SL) LA - xAZ(SL) - LB;
    xRA  = @(SL) xAZ(SL) - LA;
    ChiRA=@(SL,i) (x_i(i) <= xRA(SL)) .* exp(-(xRA(SL)-x_i(i)).^2 /Lsmooth^2) + ...
                  (x_i(i) > xRA(SL)) .* (x_i(i) < xAZ(SL)) + ...
                  (x_i(i) >= xAZ(SL)) .* exp(-(x_i(i)-xAZ(SL)).^2 /Lsmooth^2);
    ChiLA=@(SL,i) (x_i(i) <= xLA(SL)) .* exp(-(xLA(SL)-x_i(i)).^2 /Lsmooth^2) + ...
                  (x_i(i) > xLA(SL));
    Q    = @(SL) Q0 - alphaQ*(SLQ-SL)*(SL<SLQ);
    Kpn0 = Kbasic*gamma*gamma;
    Kpn1 = Kbasic*gamma*gamma;
    Knp0 = @(SL) Q(SL)*Kbasic/mu;
    Knp1 = @(SL) Q(SL)*Kbasic;
    K1on = Kon;
    K1off= Koff/mu;
    rates         = zeros(n_RU,4,4,4,4); %PC(i,a,b,c,D) = P((a,D,c)_i^t+dt | (a,b,c)_i^t ) / dt
    ones_414      = ones(n_RU,4,1,4);
    expMat        = repmat(permute([0 0 1 1; 0 0 1 1; 1 1 2 2; 1 1 2 2],[3,1,4,2]),n_RU,1,1,1);
    jMat          = repmat((1:n_RU)',1,4,1,4);
    gammaPlusExp  = gamma.^expMat;
    gammaMinusExp = gamma.^-expMat;
    rates(:,:,2,:,1) = Koff*ones_414;
    rates(:,:,3,:,4) = K1off*ones_414;
    rates(:,:,4,:,1) = Kpn0*gammaMinusExp;
    rates(:,:,3,:,2) = Kpn1*gammaMinusExp;
    P_local = @(xODE)  [sum(sum(sum(xODE(1,3:4,:,:),2),3),4); ... % first unit
                        sum(sum(sum(xODE(:,:,3:4,:),2),3),4); ... % central units
                        sum(sum(sum(xODE(end,:,:,3:4),2),3),4)];  % last unit

    %% Time loop
    times = min(input.times):dt:max(input.times);
    Ca    = interp1(input.times,input.Ca,times);
    SL    = interp1(input.times,input.SL,times);
    nT    = length(times);
    P     = zeros(1,nT);

    time_init = tic();
    fprintf('RDQ18 model. Computing... ')
    for iT = 1:nT
        if mod(iT-1, freqenceRatesUpdate) == 0
            % For performance reasons, transition rates are updated every 10 iterations
            update_rates(Ca(iT),SL(iT));
        end
        if iT == 1
            %x(i,a,b,c) = P((X_{i-1},X_i,X_{i+1}) = (a,b,c) )
            x = zeros(n_RU-2,4,4,4);
            % initial state (everything in state 1, i.e. unbounded & non-permissive)
            x(:,1,1,1) = 1;
        else
            x = x + dt * get_rhs(x);
        end
        P(iT) = mean(P_local(x));
    end
    fprintf('done. Time elapsed: %1.3f s\n', toc(time_init))

    output.times = times;
    output.Ca    = Ca;
    output.SL    = SL;
    output.P     = P;
    output.Ta    = TaMax * P;

    %% Functions definition
    function update_rates(Ca, SL)
        ChiRAmat = ChiRA(SL,jMat);
        ChiLAChiRAgammaPlus = ChiRAmat.*ChiLA(SL,jMat).*gammaPlusExp;
        rates(:,:,1,:,2) = Kon.*ChiRAmat*Ca;
        rates(:,:,2,:,3) = Knp1(SL).*ChiLAChiRAgammaPlus;
        rates(:,:,4,:,3) = K1on.*ChiRAmat*Ca;
        rates(:,:,1,:,4) = Knp0(SL).*ChiLAChiRAgammaPlus;
    end

    function rhs = get_rhs(x)
        x2_L = sum(x(2:end,:,:,:),4);                       % x2_L(i-2,a,b) = P( (X_{i-1},X_i) = (a,b) ) , i = 3:n_RU-1
        x2_R = permute(sum(x(1:end-1,:,:,:),2),[1 3 4 2]);  % x2_R(i-1,b,c) = P( (X_i,X_{i+1}) = (b,c) ) , i = 2:n_RU-2
        PhiC = rates(2:(n_RU-1),:,:,:,:) .* x; % probability fluxes (central units)
        PsiL_squeezed(2:(n_RU-2),:,:,:,:) = permute(sum(PhiC(1:(n_RU-3),:,:,:,:),2),[1,3,4,6,5,2])./ x2_L;
        PsiL_squeezed(isnan(PsiL_squeezed))=0;
        PsiL_squeezed(isinf(PsiL_squeezed))=0;
        PsiL_squeezed(1,:,:,:,:) = permute(rates(1,1,:,:,:),[1,3,4,6,5,2]);
        PhiL = PsiL_squeezed .* x; % probability fluxes (left units)
        PsiR_squeezed(1:(n_RU-3),:,:,:,:) = permute(sum(PhiC(2:(n_RU-2),:,:,:,:),4),[1,6,2,3,5,4]) ./ permute(x2_R,[1,4,2,3]);
        PsiR_squeezed(isnan(PsiR_squeezed))=0;
        PsiR_squeezed(isinf(PsiR_squeezed))=0;
        PsiR_squeezed(n_RU-2,:,:,:,:) = permute(rates(n_RU,:,:,1,:),[1,6,2,3,5,4]);
        PhiR = PsiR_squeezed .* x; % probability fluxes (right units)
        rhs = sum( permute(PhiC,[1,2,5,4,3]) - PhiC + permute(PhiL,[1,5,3,4,2]) - PhiL + permute(PhiR,[1,2,3,5,4]) - PhiR,5);
    end
end