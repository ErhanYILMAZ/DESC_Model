#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
This file is licensed under CC BY-SA (Attribution-ShareAlike)
It is provided "as is", without express or implied warranty, for educational and informational purposes only.

Python class file for Double Enhanced Self Correcting(DESC) model with 1-RC branch

Double Enhanced Self Correcting(DESC) model is an Electrical Equivalent Circuit Model 
for lithium-ion batteries also can be used for other battery types. 
Related article for DESC can be accessed here <DESC paper link will be here hopefully later>

DESC is an enhanced  derivative of Enhanced Self Correcting(ESC) Model from Gregor Plett.  
Thus, same documentation for ESC model and DESC model article can be used as referenced to understand 
code framework. After understanding the code framework and how does model work then one can develop own model.

All details, documentation and description can be found below 

http://mocha-java.uccs.edu/BMS1/index.html
http://mocha-java.uccs.edu/ECE5710/index.html
http://mocha-java.uccs.edu/ECE5710/ECE5710-Notes02.pdf
"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline



class DESC_1RC:

    def __init__(self, model, temp, dt=1.0, use_OCVS=True):

        self.temp=temp # Actual temperature
        self.model=model # Model json file for model parameters
        self.dt=dt     # Sampling period

        # if use_OCVS is false then OCV0 and OCVrel used(OCVfromSOCtempB) to calculate ocv curve
        # otherwise individual ocv curves used for test temperatures to callculate 
        # ocv curve for given temeprature (OCVfromSOCtempA)
        self.use_OCVS=use_OCVS

        # Get other parameters from the model file for the specified temperature and initialize the model
        self.Temps = np.asarray(model['Temps'])
        self.SOC = np.asarray(model['SOC'])
        self.SOC0 = np.asarray(model['SOC0'])
        self.SOCrel = np.asarray(model['SOCrel'])
        self.OCVrel = np.asarray(model['OCVrel'])
        self.OCV = np.asarray(model['OCV'])
        self.OCV0 = np.asarray(model['OCV0'])
        self.eta = self.getParamDESC('eta', temp)
        self.g = self.getParamDESC('g', temp)
        self.M = self.getParamDESC('M', temp)
        #self.OCVeta = self.getParamDESC('OCVeta', temp)
        self.Q = self.getParamDESC('Q', temp)
        self.R0cha = self.getParamDESC('R0cha', temp)
        self.R0dch = self.getParamDESC('R0dch', temp)
        self.R1 = self.getParamDESC('R1', temp)
        self.RC1 = self.getParamDESC('RC1', temp)
        self.exp1 = np.exp(-self.dt/self.RC1)
        self.one_exp1 = (1.0-self.exp1)
        self.dt_Q_3600 = self.dt/(self.Q * 3600.0) 
        self.g_dt_Q = self.g * self.dt/self.Q
        
        # Copy each OCV curve for each temperature values in the model file 
        self.OCVS= np.zeros([len(self.Temps), len(self.OCV0)])
        if use_OCVS:
            for erhan in range(len(self.Temps)):
                if self.Temps[erhan] >0:
                    self.OCVS[erhan]= model['OCV_P%02d'%self.Temps[erhan]]
                else:
                    self.OCVS[erhan]= model['OCV_N%02d'%(np.abs(self.Temps[erhan]))]
        self.ocv = self.OCVfromSOCtemp(temp)


    '''
    getParamDESC is to get model parameters for the specified temperature value
    if the temperature value is not in the tempearature values then it does inter/extrapolation
    to obtain parameters.
    '''
    def getParamDESC(self, paramName,temp):
        if not np.isscalar(temp):
            raise Exception('Temp has to be scalar!')
    
        if temp in self.Temps:
            index, =np.where(self.Temps == temp)
            return self.model[paramName][index[0]]
        else:
            # do inter/extrapolation
            # k spline order: 1 linear, 2 quadratic, 3 cubic ... 
            s = InterpolatedUnivariateSpline(self.Temps, self.model[paramName], k=1)
            return s(temp)


    '''
    parameter_update additionally can update model parameters later
    the model instant created for specified temperature value.
    '''
    def parameter_update(self, temp):
        self.eta = self.getParamDESC('eta', temp)
        self.g = self.getParamDESC('g', temp)
        self.M = self.getParamDESC('M', temp)
        #self.OCVeta = self.getParamDESC('OCVeta', temp)
        self.Q = self.getParamDESC('Q', temp)
        self.R0cha = self.getParamDESC('R0cha', temp)
        self.R0dch = self.getParamDESC('R0dch', temp)
        self.R1 = self.getParamDESC('R1', temp)
        self.RC1 = self.getParamDESC('RC1', temp)
        self.ocv = self.OCVfromSOCtemp(temp)
        self.temp=temp
        self.exp1 = np.exp(-self.dt/self.RC1)
        self.one_exp1 = (1.0-self.exp1)
        self.dt_Q_3600 = self.dt/(self.Q * 3600.0) 
        self.g_dt_Q = self.g * self.dt/self.Q


    '''
    state_equation is state transition equation for DESC model
    It calculates next state for given previous state and current as input
    It convenient for recursive filtering
    '''
    def state_equation(self,state,dt, current):
        
        
        # Calculate exponential term for RC branch state
        
        #if current < 0.0:
        #    current = current * self.eta
        
        current = np.where(current<=0, current * self.eta, current )

        # Calculate exponential term for hysteresis state
        exp2 = np.exp(-(np.abs(current) * self.g_dt_Q))

        # SoC state calculation for next time step
        state[0] = state[0] - self.dt_Q_3600* current  

        # RC branch state calculation for next time step
        state[1] = self.exp1 * state[1] + self.one_exp1 * current

        # Hysteresis state calculation for next time step
        state[2] = exp2 * state[2] - (exp2-1.0) * np.sign(current)
        # Return calculated states
        
        # Return calculated next state
        return state


    '''
    output_equation is an output function for DESC model
    It calculates output of DESC model for given state and current as input
    It convenient for recursive filtering
    '''
    def output_equation(self,state, current):
        
        # Calculate vr0 value according to imput current polarity
        #if current < 0.0:
        #    current = current * self.eta
        #    vr0 = self.R0dch * current
        #else:
        #    vr0 = self.R0cha  * current
            
        # The line bbelow does the same as above. It can be faster for batch calculation
        vr0 = np.where(current<=0,self.R0dch * current * self.eta ,self.R0cha * current )

        # Calculate OCV value for the given soc value
        OCV = self.OCVfromSOC(state[0], self.ocv)
 
        # Calculate cell output voltage
        vout = OCV + self.M * state[2] - self.R1 * state[1] - vr0

        #return output voltage
        return vout


    '''
    It calculates directly output of DESC model for given previous state and current as input
    It can be used for batch runnin or simulation
    '''
    def fun(self, state, current):
        vout=np.zeros(current.size)

        for i in range(current.size):
            state = self.state_equation(state,self.dt,current[i])
            vout[i]=self.output_equation(state,current[i])
        return vout


    '''
    SOCfromOCV is works as look up table and it calculates soc value from the ocv
    it is implemented from here https://github.com/batterysim/esctoolbox-python
    '''
    def SOCfromOCV(self, ocv):
        """ OCV function """
        OCV = np.asarray(self.ocv)        # force to be column vector  
        SOC = np.asarray(self.SOC)          # force to be column vector
        ocvcol = np.asarray(ocv)
        # do inter/extrapolation
        
        if ocvcol.ndim == 0:
            ocvcol = ocvcol[None]
        
        diffOCV = OCV[1] - OCV[0]           # spacing between SOC points - assume uniform
        soc = np.zeros(np.size(ocvcol))     # initialize output to zero
        I1, = np.where(ocvcol <= OCV[0])    # indices of socs below model-stored data
        I2, = np.where(ocvcol >= OCV[-1])   # and of socs above model-stored data
        I3, = np.where((ocvcol > OCV[0]) & (ocvcol < OCV[-1]))   # the rest of them
        I6 = np.isnan(ocvcol)               # if input is "not a number" for any locations

    # for voltages less than lowest stored soc datapoint, extrapolate off
    # low end of table    
        if I1.size:
            soc[I1] = SOC[0] - (SOC[1]-OCV[0])*(OCV[0]-ocvcol[I1])/(OCV[1]-OCV[0])
            
    # for voltages greater than highest stored soc datapoint, extrapolate off
    # high end of table        
        if I2.size:
            soc[I2] = SOC[-1] + (ocvcol[I2]-OCV[-1])*(SOC[-1]-SOC[-2])/(OCV[-1]-OCV[-2])
            
        # for normal soc range, manually interpolate (10x faster than "interp1")
        I4 = (ocvcol[I3] - OCV[0])/diffOCV  # using linear interpolation
        I5 = np.floor(I4)
        I5 = I5.astype(int)
        I45 = I4 - I5
        omI45 = 1 - I45
        soc[I3] = SOC[I5]*omI45 + SOC[I5+1]*I45
        soc[I6] = 0     # replace NaN SOCs with zero voltage
        return soc


    '''
    OCVfromSOC is works as look up table and it calculates ocv value from the soc
    it is implemented from here https://github.com/batterysim/esctoolbox-python
    '''
    def OCVfromSOC(self,soc, OCV):
    
        """ OCV function """
        OCV = np.asarray(OCV)        # force to be column vector  
        SOC = np.asarray(self.SOC)          # force to be column vector
        soccol = np.asarray(soc)
        
        
        if soccol.ndim == 0:
            soccol = soccol[None]
        
        diffSOC = SOC[1] - SOC[0]           # spacing between SOC points - assume uniform
        ocv = np.zeros(np.size(soccol))     # initialize output to zero
        I1, = np.where(soccol <= SOC[0])    # indices of socs below model-stored data
        I2, = np.where(soccol >= SOC[-1])   # and of socs above model-stored data
        I3, = np.where((soccol > SOC[0]) & (soccol < SOC[-1]))   # the rest of them
        I6 = np.isnan(soccol)               # if input is "not a number" for any locations

    # for voltages less than lowest stored soc datapoint, extrapolate off
    # low end of table    
        if I1.size:
            ocv[I1] = OCV[0] - (OCV[1]-OCV[0])*(SOC[0]-soccol[I1])/(SOC[1]-SOC[0])
            
    # for voltages greater than highest stored soc datapoint, extrapolate off
    # high end of table        
        if I2.size:
            ocv[I2] = OCV[-1] + (soccol[I2]-SOC[-1])*(OCV[-1]-OCV[-2])/(SOC[-1]-SOC[-2])
            
        # for normal soc range, manually interpolate (10x faster than "interp1")
        I4 = (soccol[I3] - SOC[0])/diffSOC  # using linear interpolation
        I5 = np.floor(I4)
        I5 = I5.astype(int)
        I45 = I4 - I5
        omI45 = 1 - I45
        ocv[I3] = OCV[I5]*omI45 + OCV[I5+1]*I45
        ocv[I6] = 0     # replace NaN SOCs with zero voltage
        return ocv


    '''
    This function calculates OCV curves for the desired temperatures from a given SoC values.
    '''
    def OCVfromSOCtemp(self, temp):
        if self.use_OCVS == True:
            return self.OCVfromSOCtempA(temp)
        else:
            return self.OCVfromSOCtempB(temp)


    '''
    OCVfromSOCtempA uses OCV values those obtained from test
    temperature values to calculate OCV value for the specified temperature value.
    If the temperature is in the test temparature values then it returns directly the
    related OCV curve otherwise it uses interpolation/extrapolation to calculate
    OCV values for intermediate, lower or higher temperature values.
    '''
    def OCVfromSOCtempA(self, temp):
        if temp in self.Temps:
            index, =np.where(self.Temps == temp)
            return self.OCVS[index[0]]
        else:
            # do inter/extrapolation
            # k spline order: 1 linear, 2 quadratic, 3 cubic ... 
            s=[0]*self.OCVS.shape[1]
            for erhan in range(self.OCVS.shape[1]):
                e = InterpolatedUnivariateSpline(self.Temps, self.OCVS[:,erhan], k=3)
                s[erhan] =e(temp)
            return s


    '''
    This function uses OCV0 and OCVrel to calculate OCV for the specified temperature value
    As I experienced to get OCV0 and OCVrel from OCV values then compute OCV values lower 
    the accuracy. Thus, I have deprecated this method and instead used OCVfromSOCtemp to calculate
    OCV values for the specified temperature values.
    '''
    def OCVfromSOCtempB(self, temp):
        OCV0 = self.OCVfromSOC(self.SOC, self.OCV0)
        OCVrel = self.OCVfromSOC(self.SOC, self.OCVrel)
        return OCV0 + temp * OCVrel
