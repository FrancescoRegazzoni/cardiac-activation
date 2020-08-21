#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

import numpy as np
from scipy.linalg import expm
import json
import time

class model_RDQ20_SE:
    """
    Class implementing the spatially-explicitly ODE model (SE-ODE) for
    cardiomyocytes force generation presented in [1, 2].

    References
    ----------

    [1] F. Regazzoni "Mathematical modeling and Machine Learning for the
        numerical simulation of cardiac electromechanics", PhD Thesis -
        Politecnico di Milano (2020)
        http://hdl.handle.net/10589/152617
    [2] F. Regazzoni, L. Dede', A. Quarteroni "Biophysically detailed
        mathematical models of multiscale cardiac active mechanics",
        submitted (2020)
        https://arxiv.org/abs/2004.07910
    """

    def __init__(self, params = '../params/params_RDQ20-SE_human_body-temperature.json'):
        """
        Constructor.

        Parameters
        ----------
        params : dictionary containing the model parameters, or path of a json
                 file
        """

        if isinstance(params, str):
            with open(params) as json_file:
                params = json.load(json_file)

        self.LA       = float(params['geometry']['LA'])             # [micro m]
        self.LM       = float(params['geometry']['LM'])             # [micro m]
        self.LB       = float(params['geometry']['LB'])             # [micro m]
        self.SL0      = float(params['geometry']['SL0'])            # [micro m]
        self.n_RU     =   int(params['geometry']['n_RU'])           # [-]
        self.n_MH     =   int(params['geometry']['n_MH'])           # [-]
        self.Lsmooth  = float(params['geometry']['Lsmooth'])        # [micro m]
        self.Q        = float(params['RU_steady_state']['Q'])       # [-]
        self.Kd0      = float(params['RU_steady_state']['Kd0'])     # [micro M]
        self.alphaKd  = float(params['RU_steady_state']['alphaKd']) # [micro M / micro m]
        self.mu       = float(params['RU_steady_state']['mu'])      # [-]
        self.gamma    = float(params['RU_steady_state']['gamma'])   # [-]
        self.Koff     = float(params['RU_kinetics']['Koff'])        # [s^-1]
        self.Kbasic   = float(params['RU_kinetics']['Kbasic'])      # [s^-1]
        self.r0       = float(params['XB_cycling']['r0'])           # [s^-1]
        self.alpha    = float(params['XB_cycling']['alpha'])        # [-]
        self.mu0_fP   = float(params['XB_cycling']['mu0_fP'])       # [s^-1]
        self.mu1_fP   = float(params['XB_cycling']['mu1_fP'])       # [s^-1]
        self.a_XB     = float(params['upscaling']['a_XB'])          # [kPa]

        # Numerical parameters
        self.dt_RU = 2.5e-5                                     # [s]
        self.freqXB = round(1e-3 / self.dt_RU)                  # [-]
        self.freqenceRatesUpdate = round(2.5e-4 / 2.5e-5)       # [-]

        # Initialization
        self.y_j = lambda j: self.LA * (j-.5)/self.n_RU
        self.yLA = lambda SL: 2*self.LA - SL
        self.yM0 = lambda SL: (2*self.LA - SL + self.LB) / 2
        self.yM1 = lambda SL: (2*self.LA - SL + self.LM) / 2
        self.ChiMF = lambda SL, j: .5*np.tanh( (self.y_j(j) - self.yM0(SL))/self.Lsmooth) \
                                 + .5*np.tanh(-(self.y_j(j) - self.yM1(SL))/self.Lsmooth)
        self.ChiSF = lambda SL, j: .5*np.tanh((self.y_j(j) - self.yLA(SL))/self.Lsmooth) + .5
        self.Kd = lambda SL: self.Kd0 - self.alphaKd*(2.15 - SL)

        self.Kpn0 = self.Kbasic * self.gamma * self.gamma
        self.Kpn1 = self.Kbasic * self.gamma * self.gamma
        self.Knp0 = self.Q * self.Kbasic / self.mu
        self.Knp1 = self.Q * self.Kbasic
        self.K1off= self.Koff / self.mu

        self.ratesRU = np.zeros((self.n_RU,4,4,4,4)) # rates(i,a,b,c,D) = P((a,D,c)_i^t+dt | (a,b,c)_i^t ) / dt
        self.PhiL = np.zeros((self.n_RU-2,4,4,4,4)) # (i,a,b,c,D)
        self.PhiR = np.zeros((self.n_RU-2,4,4,4,4)) # (i,a,b,c,D)
        self.jMat = np.arange(1,self.n_RU+1)
        self.expMat = np.array([[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 2, 2],[1, 1, 2, 2]])[None,:,:]
        self.gammaPlusExp = self.gamma**self.expMat
        self.gammaMinusExp = self.gamma**-self.expMat
        self.ratesRU[:,:,1,:,0] = self.Koff
        self.ratesRU[:,:,3,:,0] = self.Kpn0*self.gammaMinusExp
        self.ratesRU[:,:,2,:,1] = self.Kpn1*self.gammaMinusExp
        self.ratesRU[:,:,2,:,3] = self.K1off
        self.P_local = lambda x: np.concatenate((x[0,2:4,:,:][None,:].sum(axis=(1,2,3)),    # first unit
                                                 x[:,:,2:4,:].sum(axis=(1,2,3)),            # central units
                                                 x[-1,:,:,2:4][None,:].sum(axis=(1,2,3))))  # last unit

    def solve(self, inputs):
        """
        Perform a simulation with the model.

        Parameters
        ----------
        inputs: dictionary containing the input data
          - inputs['times']: time instants [s]
          - inputs['Ca']:    intracellular calcium ions concentration [micro M]
          - inputs['SL']:    sarcomere length [micro m]

        Returns
        -------
        output: dictionary containing the output data
          - output['times']: time instants [s]
          - output['Ca']:    intracellular calcium ions concentration [micro M]
          - output['SL']:    sarcomere length [micro m]
          - output['P']:     permissivity [-]
          - output['Ta']:    active tension [kPa]
        """

        times = np.arange(np.min(inputs['times']), np.max(inputs['times']), self.dt_RU)
        Ca    = np.interp(times, inputs['times'], inputs['Ca'])
        SL    = np.interp(times, inputs['times'], inputs['SL'])
        nT    = len(times)
        P     = np.zeros(nT)
        Ta    = np.zeros(nT)

        old_settings = np.seterr(invalid='ignore') # ignore divide-by-zero warning
        time_init = time.time();
        print('RDQ20-SE model. Computing... ', end = '')
        for iT in range(nT):
            if iT % self.freqenceRatesUpdate == 0:
                # For performance reasons, transition rates are updated every 10 iterations
                self.__RU_update_rates(Ca[iT],SL[iT])
            if iT == 0:
                # x_RU(i,a,b,c) = P((X_{i-1},X_i,X_{i+1}) = (a,b,c) )
                x_RU = np.zeros((self.n_RU-2,4,4,4))
                # initial state (everything in state 1, i.e. unbounded & non-permissive)
                x_RU[:,0,0,0] = 1.
                x_XB = np.zeros((self.n_MH,2,2))
            else:
                x_RU = x_RU + self.dt_RU * self.__RU_get_rhs(x_RU)
                if iT % self.freqXB == 0:
                    dSLdt = (SL[iT] - SL[iT - self.freqXB]) / (self.dt_RU * self.freqXB)
                    x_XB = self.__XB_advance(x_RU, x_XB, self.dt_RU * self.freqXB, SL[iT], dSLdt)
            P[iT]  = np.mean(self.P_local(x_RU))
            Ta[iT] = self.a_XB * np.mean(np.sum(x_XB[:,1,:], axis = 1))
        print('done. Time elapsed: %1.3f s' % (time.time() - time_init))
        np.seterr(**old_settings)

        output = {}
        output['times'] = times
        output['Ca']    = Ca
        output['SL']    = SL
        output['P']     = P
        output['Ta']    = Ta

        return output

    def __RU_update_rates(self, Ca, SL):
        Kon = self.Koff / self.Kd(SL) # [micro M^-1 * s^-1]
        Chi_2_gammaPlus = self.ChiSF(SL,self.jMat)[:,None,None]*self.gammaPlusExp
        self.ratesRU[:,:,0,:,1] = Kon*Ca
        self.ratesRU[:,:,1,:,2] = self.Knp1*Chi_2_gammaPlus
        self.ratesRU[:,:,3,:,2] = Kon*Ca
        self.ratesRU[:,:,0,:,3] = self.Knp0*Chi_2_gammaPlus

    def __RU_get_rhs(self, x_RU):
        x2 = x_RU.sum(axis=3)
        self.PhiC = self.ratesRU[1:-1,:,:,:,:] * x_RU[:,:,:,:,None] # probability fluxes (central units)
        self.PhiL[1:,:,:,:,:] = self.PhiC[:-1,:,:,:,:].sum(axis=1)[:,:,:,None,:] / x2[1:,:,:,None,None]
        self.PhiL[np.logical_or(np.isnan(self.PhiL), np.isinf(self.PhiL))] = .0
        self.PhiL[0,:,:,:,:] = self.ratesRU[0,0,:,:,None,:]
        self.PhiL = self.PhiL * x_RU[:,:,:,:,None] # probability fluxes (left units)
        self.PhiR[:-1,:,:,:,:] = self.PhiC[1:,:,:,:,:].sum(axis=3)[:,None,:,:,:] / x2[1:,None,:,:,None]
        self.PhiR[np.logical_or(np.isnan(self.PhiR), np.isinf(self.PhiR))] = .0
        self.PhiR[-1,:,:,:,:] = self.ratesRU[-1,None,:,:,0,:]
        self.PhiR = self.PhiR * x_RU[:,:,:,:,None] # probability fluxes (right units)
        return ( self.PhiC.swapaxes(2,4) - self.PhiC \
               + self.PhiL.swapaxes(1,4) - self.PhiL \
               + self.PhiR.swapaxes(3,4) - self.PhiR ).sum(axis=4)

    def __XB_advance(self, x_RU, x_XB, dt, SL, dSLdt):
        v = -dSLdt / self.SL0 # [1/s]
        x_RU_new = np.zeros((self.n_MH,2,2))

        perm = self.P_local(x_RU)
        for i in range(2):
            if i == 0:
                a1, a2, b1, b2 = 2, 3, 1, 0
                a_marginal = perm
            else:
                a1, a2, b1, b2 = 1, 0, 2, 3
                a_marginal = 1.0 - perm
            k_L = np.sum(self.ratesRU[0,0,a1,:,b1][:,None] * x_RU[0,a1,:,:], axis = (0,1)) + \
                  np.sum(self.ratesRU[0,0,a2,:,b2][:,None] * x_RU[0,a2,:,:], axis = (0,1))
            k_C = np.sum(self.ratesRU[1:-1,:,a1,:,b1]*x_RU[:,:,a1,:], axis = (1,2)) + \
                  np.sum(self.ratesRU[1:-1,:,a2,:,b2]*x_RU[:,:,a2,:], axis = (1,2))
            k_R = np.sum(self.ratesRU[-1,:,a1,0,b1]*x_RU[-1,:,:,a1], axis = (0,1)) + \
                  np.sum(self.ratesRU[-1,:,a2,0,b2]*x_RU[-1,:,:,a2], axis = (0,1))
            k_perm = np.concatenate((np.array(k_L)[None], k_C, np.array(k_R)[None])) / a_marginal
            if i == 0:
                k_PN = k_perm
            else:
                k_NP = k_perm

        chiSF_i = self.ChiSF(SL, np.arange(1,self.n_RU+1)) * self.ChiMF(SL, np.arange(1,self.n_RU+1))
        mu0_fP_i = self.mu0_fP*chiSF_i
        mu1_fP_i = self.mu1_fP*chiSF_i
        r = self.r0 + self.alpha*abs(v)
        diag_P = r + k_PN
        diag_N = r + k_NP

        for i in range(self.n_MH):
            sol_i = np.reshape(x_XB[i,:,:], (4,), order='F')
            rhs = np.array((perm[i]*mu0_fP_i[i], perm[i]*mu1_fP_i[i], 0.0, 0.0))
            A = np.array([[-diag_P[i], 0         , k_NP[i]   , 0         ], \
                          [-v        , -diag_P[i], 0         , k_NP[i]   ], \
                          [k_PN[i]   , 0         , -diag_N[i], 0         ], \
                          [0         , k_PN[i]   , -v        , -diag_N[i]]])
            sol_inf_i = -np.linalg.solve(A, rhs)
            x_RU_new[i,:,:] = np.reshape(sol_inf_i + np.matmul(expm(dt*A), sol_i - sol_inf_i), (2,2), order='F')

        return x_RU_new