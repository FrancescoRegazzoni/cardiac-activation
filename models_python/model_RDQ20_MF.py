#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

import numpy as np
from scipy.linalg import expm
import json
import time

class model_RDQ20_MF:
    """
    Class implementing the mean-field ODE model (MF-ODE) for cardiomyocytes
    force generation presented in [1, 2].

    References
    ----------

    [1] F. Regazzoni "Mathematical modeling and Machine Learning for the
        numerical simulation of cardiac electromechanics", PhD Thesis -
        Politecnico di Milano (2020)
        http://hdl.handle.net/10589/152617
    [2] F. Regazzoni, L. Dede', A. Quarteroni "Biophysically detailed
        mathematical models of multiscale cardiac active mechanics",
        PLOS Computational Biology (2020)
        https://doi.org/10.1371/journal.pcbi.1008294
    """

    def __init__(self, params = '../params/params_RDQ20-MF_human_body-temperature.json'):
        """
        Constructor.

        Parameters
        ----------
        params : dictionary containing the model parameters, or path of a json
                 file
        """

        self.model_name = 'RDQ20-MF'

        if isinstance(params, str):
            with open(params) as json_file:
                params = json.load(json_file)

        self.LA       = float(params['geometry']['LA'])             # [micro m]
        self.LM       = float(params['geometry']['LM'])             # [micro m]
        self.LB       = float(params['geometry']['LB'])             # [micro m]
        self.SL0      = float(params['geometry']['SL0'])            # [micro m]
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
        LMh = (self.LM - self.LB) * .5;
        self.frac_SO =  lambda SL: ((SL>self.LA)          *(SL<=self.LM)          *(SL-self.LA) + \
                                    (SL>self.LM)          *(SL<=2*self.LA-self.LB)*(SL+self.LM-2*self.LA)*.5 + \
                                    (SL>2*self.LA-self.LB)*(SL<=2*self.LA+self.LB)*LMh + \
                                    (SL>2*self.LA+self.LB)*(SL<=2*self.LA+self.LM)*(self.LM+2*self.LA-SL)*.5 \
                                   )/LMh
        self.Kd = lambda SL: self.Kd0 - self.alphaKd*(2.15 - SL)
        self.kC = np.zeros((2,2))     # kC(A,a) = P( C_i^t+dt =~ A | (C_i,T_i)^t = (A,a) ) / dt
        self.kT = np.zeros((2,2,2,2)) # kT(a,b,c,B) = P( T_i^t+dt =~ b | (T_i-1,T_i,T_i+1,C_i)^t = (a,b,c,B) ) / dt
        expMat = np.array([[0, 1],[1, 2]])
        self.kC[1,0] = self.Koff
        self.kC[1,1] = self.Koff / self.mu
        self.kT[:,1,:,0] = self.Kbasic * self.gamma**(2 - expMat)
        self.kT[:,1,:,1] = self.Kbasic * self.gamma**(2 - expMat)
        self.kT[:,0,:,0] = self.Q * self.Kbasic / self.mu * self.gamma**expMat
        self.kT[:,0,:,1] = self.Q * self.Kbasic * self.gamma**expMat

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
        print('RDQ20-MF model. Computing... ', end = '')
        for iT in range(nT):
            if iT % self.freqenceRatesUpdate == 0:
                # For performance reasons, transition rates are updated every 10 iterations
                self.__RU_update_rates(Ca[iT],SL[iT])
            if iT == 0:
                # x_RU(a,b,c,B) = P( (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) )
                x_RU = np.zeros((2,2,2,2))
                # initial state (everything in state 1, i.e. unbounded & non-permissive)
                x_RU[0,0,0,0] = 1.
                x_XB = np.zeros((2,2))
            else:
                x_RU = x_RU + self.dt_RU * self.__RU_get_rhs(x_RU)
                if iT % self.freqXB == 0:
                    dSLdt = (SL[iT] - SL[iT - self.freqXB]) / (self.dt_RU * self.freqXB)
                    x_XB = self.__XB_advance(x_RU, x_XB, self.dt_RU * self.freqXB, dSLdt)
            P[iT]  = x_RU[:,1,:,:].sum() * self.frac_SO(SL[iT])
            Ta[iT] = self.a_XB * np.sum(x_XB[1,:]) * self.frac_SO(SL[iT])
        print('done. Time elapsed: %1.3f s' % (time.time() - time_init))
        np.seterr(**old_settings)

        output = {}
        output['t']  = times
        output['Ca'] = Ca
        output['SL'] = SL
        output['P']  = P
        output['Ta'] = Ta

        return output

    def __RU_update_rates(self, Ca, SL):
        Kon = self.Koff / self.Kd(SL) # [micro M^-1 * s^-1]
        self.kC[0,:] = Kon*Ca

    def __RU_get_rhs(self, x_RU):
        PhiT_C = np.nan_to_num(x_RU*self.kT, nan = 0.0)                                # ...(a,b,c,B) =  P( T_i^t+dt =~ b, (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        kT_L   = np.sum(PhiT_C, axis = (0,3)) / np.sum(x_RU, axis = (0,3))             # ...(a,b) ~= P( T_{i-1}^t+dt =~ a | (T_{i-1},T_i)^t = (a,b) ) / dt
        kT_R   = np.sum(PhiT_C, axis = (2,3)) / np.sum(x_RU, axis = (2,3))             # ...(b,c) ~= P( T_{i+1}^t+dt =~ c | (T_i,T_{i+1})^t = (b,c) ) / dt
        PhiT_L = np.nan_to_num(kT_L[:,:,None,None]*x_RU, nan = 0.0)                    # ...(a,b,c,B) ~= P( T_{i-1}^t+dt =~ a , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        PhiT_R = np.nan_to_num(kT_R[None,:,:,None]*x_RU, nan = 0.0)                    # ...(a,b,c,B) ~= P( T_{i+1}^t+dt =~ c , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        PhiC_C = np.nan_to_num(x_RU * self.kC.swapaxes(1,0)[None,:,None,:], nan = 0.0) # ...(a,b,c,B) =  P( C_i^t+dt =~ B     , (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
        return - PhiT_L + np.flip(PhiT_L, axis = 0) \
               - PhiT_C + np.flip(PhiT_C, axis = 1) \
               - PhiT_R + np.flip(PhiT_R, axis = 2) \
               - PhiC_C + np.flip(PhiC_C, axis = 3)

    def __XB_advance(self, x_RU, x_XB, dt, dSLdt):
        v = -dSLdt / self.SL0 # [1/s]
        perm = x_RU[:,1,:,:].sum()
        k_PN = np.nan_to_num((self.kT[:,1,:,:] * x_RU[:,1,:,:]).sum() / x_RU[:,1,:,:].sum(), nan = 0.0)
        k_NP = np.nan_to_num((self.kT[:,0,:,:] * x_RU[:,0,:,:]).sum() / x_RU[:,0,:,:].sum(), nan = 0.0)
        r = self.r0 + self.alpha*abs(v)
        diag_P = r + k_PN
        diag_N = r + k_NP
        sol = np.reshape(x_XB, (4,), order='F')
        rhs = np.array((perm*self.mu0_fP, perm*self.mu1_fP, 0.0, 0.0))
        A = np.array([[-diag_P, 0      , k_NP   , 0      ], \
                      [-v     , -diag_P, 0      , k_NP   ], \
                      [k_PN   , 0      , -diag_N, 0      ], \
                      [0      , k_PN   , -v     , -diag_N]])
        sol_inf = -np.linalg.solve(A, rhs)
        return np.reshape(sol_inf + np.matmul(expm(dt*A), sol - sol_inf), (2,2), order='F')