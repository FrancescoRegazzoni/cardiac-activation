#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

import numpy as np
import json
import time

class model_RDQ18:
    """
    Class implementing the ODE model for sarcomere activation proposed in [1].

    References
    ----------

    [1] F. Regazzoni, L. Dedè, A. Quarteroni "Active contraction of cardiac
        cells: a reduced model for sarcomere dynamics with cooperative
        interactions", Biomechanics and Modeling in Mechanobiology (2018)
        https://doi.org/10.1007/s10237-018-1049-0
    """

    def __init__(self, params = '../params/params_RDQ18.json'):
        """
        Constructor.

        Parameters
        ----------
        params : dictionary containing the model parameters, or path of a json
                 file
        """

        self.model_name = 'RDQ18'

        if isinstance(params, str):
            with open(params) as json_file:
                params = json.load(json_file)

        self.LA       = float(params['geometry']['LA'])             # [micro m]
        self.LM       = float(params['geometry']['LM'])             # [micro m]
        self.LB       = float(params['geometry']['LB'])             # [micro m]
        self.n_RU     =   int(params['geometry']['n_RU'])           # [-]
        self.Lsmooth  = float(params['geometry']['Lsmooth'])        # [micro m]
        self.Q0       = float(params['RU_steady_state']['Q0'])      # [-]
        self.SLQ      = float(params['RU_steady_state']['SLQ'])     # [micro m]
        self.alphaQ   = float(params['RU_steady_state']['alphaQ'])  # [micro m^-1]
        self.mu       = float(params['RU_steady_state']['mu'])      # [-]
        self.gamma    = float(params['RU_steady_state']['gamma'])   # [-]
        self.Kon      = float(params['RU_kinetics']['Kon'])         # [micro M^-1 * s^-1]
        self.Koff     = float(params['RU_kinetics']['Koff'])        # [s^-1]
        self.Kbasic   = float(params['RU_kinetics']['Kbasic'])      # [s^-1]
        self.TaMax    = float(params['upscaling']['TaMax'])         # [kPa]

        # Numerical parameters
        self.dt = 2.5e-5                         # [s]
        self.freqenceRatesUpdate = 10            # [-]

        # Initialization
        self.x_i = lambda i: (self.LM-self.LB)*.5 * i/self.n_RU
        self.xAZ = lambda SL: (SL-self.LB)/2.
        self.xLA = lambda SL: self.LA -self.xAZ(SL) - self.LB
        self.xRA = lambda SL: self.xAZ(SL) - self.LA
        self.ChiRA = lambda SL,i: (self.x_i(i) <= self.xRA(SL)) * np.exp(-(self.xRA(SL)-self.x_i(i))**2 /self.Lsmooth**2) + \
                      (self.x_i(i) > self.xRA(SL)) * (self.x_i(i) < self.xAZ(SL)) + \
                      (self.x_i(i) >= self.xAZ(SL)) * np.exp(-(self.x_i(i)-self.xAZ(SL))**2 /self.Lsmooth**2)
        self.ChiLA = lambda SL,i: (self.x_i(i) <= self.xLA(SL)) * np.exp(-(self.xLA(SL)-self.x_i(i))**2 /self.Lsmooth**2) + \
                      (self.x_i(i) > self.xLA(SL))
        self.Q = lambda SL: self.Q0 - self.alphaQ*(self.SLQ-SL)*(SL<self.SLQ)
        self.Kpn0 = self.Kbasic*self.gamma*self.gamma
        self.Kpn1 = self.Kbasic*self.gamma*self.gamma
        self.Knp0 = lambda SL: self.Q(SL)*self.Kbasic/self.mu
        self.Knp1 = lambda SL: self.Q(SL)*self.Kbasic
        self.K1on = self.Kon
        self.K1off = self.Koff/self.mu

        self.rates = np.zeros((self.n_RU,4,4,4,4)) # rates(i,a,b,c,D) = P((a,D,c)_i^t+dt | (a,b,c)_i^t ) / dt
        self.PhiL = np.zeros((self.n_RU-2,4,4,4,4)) # (i,a,b,c,D)
        self.PhiR = np.zeros((self.n_RU-2,4,4,4,4)) # (i,a,b,c,D)
        self.jMat = np.arange(1,self.n_RU+1)
        self.expMat = np.array([[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 2, 2],[1, 1, 2, 2]])[None,:,:]
        self.gammaPlusExp = self.gamma**self.expMat
        self.gammaMinusExp = self.gamma**-self.expMat
        self.rates[:,:,1,:,0] = self.Koff
        self.rates[:,:,3,:,0] = self.Kpn0*self.gammaMinusExp
        self.rates[:,:,2,:,1] = self.Kpn1*self.gammaMinusExp
        self.rates[:,:,2,:,3] = self.K1off
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

        times = np.arange(np.min(inputs['times']), np.max(inputs['times']), self.dt)
        Ca    = np.interp(times, inputs['times'], inputs['Ca'])
        SL    = np.interp(times, inputs['times'], inputs['SL'])
        nT    = len(times)
        P     = np.zeros(nT)

        old_settings = np.seterr(invalid='ignore') # ignore divide-by-zero warning
        time_init = time.time();
        print('RDQ18 model. Computing... ', end = '')
        for iT in range(nT):
            if iT % self.freqenceRatesUpdate == 0:
                # For performance reasons, transition rates are updated every 10 iterations
                self.__update_rates(Ca[iT],SL[iT])
            if iT == 0:
                #x(i,a,b,c) = P((X_{i-1},X_i,X_{i+1}) = (a,b,c) )
                x = np.zeros((self.n_RU-2,4,4,4))
                # initial state (everything in state 1, i.e. unbounded & non-permissive)
                x[:,0,0,0] = 1.
            else:
                x = x + self.dt * self.__get_rhs(x)
            P[iT] = np.mean(self.P_local(x))
        print('done. Time elapsed: %1.3f s' % (time.time() - time_init))
        np.seterr(**old_settings)

        output = {}
        output['t']  = times
        output['Ca'] = Ca
        output['SL'] = SL
        output['P']  = P
        output['Ta'] = self.TaMax * P

        return output

    def __update_rates(self, Ca, SL):
        ChiRAmat = self.ChiRA(SL,self.jMat)
        ChiLAChiRAgammaPlus = ChiRAmat[:,None,None]*self.ChiLA(SL,self.jMat)[:,None,None]*self.gammaPlusExp
        self.rates[:,:,0,:,1] = self.Kon*ChiRAmat[:,None,None]*Ca
        self.rates[:,:,1,:,2] = self.Knp1(SL)*ChiLAChiRAgammaPlus
        self.rates[:,:,3,:,2] = self.K1on*ChiRAmat[:,None,None]*Ca
        self.rates[:,:,0,:,3] = self.Knp0(SL)*ChiLAChiRAgammaPlus

    def __get_rhs(self, x):
        x2 = x.sum(axis=3)
        self.PhiC = self.rates[1:-1,:,:,:,:] * x[:,:,:,:,None] # probability fluxes (central units)
        self.PhiL[1:,:,:,:,:] = self.PhiC[:-1,:,:,:,:].sum(axis=1)[:,:,:,None,:] / x2[1:,:,:,None,None]
        self.PhiL[np.logical_or(np.isnan(self.PhiL), np.isinf(self.PhiL))] = .0
        self.PhiL[0,:,:,:,:] = self.rates[0,0,:,:,None,:]
        self.PhiL = self.PhiL * x[:,:,:,:,None] # probability fluxes (left units)
        self.PhiR[:-1,:,:,:,:] = self.PhiC[1:,:,:,:,:].sum(axis=3)[:,None,:,:,:] / x2[1:,None,:,:,None]
        self.PhiR[np.logical_or(np.isnan(self.PhiR), np.isinf(self.PhiR))] = .0
        self.PhiR[-1,:,:,:,:] = self.rates[-1,None,:,:,0,:]
        self.PhiR = self.PhiR * x[:,:,:,:,None] # probability fluxes (right units)
        return ( self.PhiC.swapaxes(2,4) - self.PhiC \
               + self.PhiL.swapaxes(1,4) - self.PhiL \
               + self.PhiR.swapaxes(3,4) - self.PhiR ).sum(axis=4)