#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
This file is licensed under CC BY-SA (Attribution-ShareAlike)
It is provided "as is", without express or implied warranty, for educational and informational purposes only.

Python file for Double Enhanced Self Correcting(DESC) model with 1-RC branch parameter estimation.
This file does paramater estimation for model parameters from dynamic test data and uses nonlinear least squares
optimisation for parameter estimation.

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


import os
import time
import json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from DESC_1RC import DESC_1RC
from datetime import datetime

now = datetime.now()
start = now.strftime("%d/%m/%Y %H:%M:%S")

temps=[-5, 5, 15, 25, 35, 45]
mags= [20, 20, 40, 50, 50, 50]

data=[0]*len(temps)


# Read UDDS dynamic data for all temperature values
for erhan in range(len(temps)):
    
    data[erhan] = {'script1':0, 'script2':0, 'script3':0}
    if(temps[erhan]>0):
        script1  = "dyn_data/THUN100_DYN_%02d_P%02d_S1.csv" %(mags[erhan], temps[erhan]) 
        data[erhan]['script1'] = pd.read_csv(script1)
        
        script2 = "dyn_data/THUN100_DYN_%02d_P%02d_S2.csv" %(mags[erhan], temps[erhan])
        data[erhan]['script2'] = pd.read_csv(script2)
        
        script3 = "dyn_data/THUN100_DYN_%02d_P%02d_S3.csv" %(mags[erhan], temps[erhan])  
        data[erhan]['script3'] = pd.read_csv(script3)
    else:
        script1  = "dyn_data/THUN100_DYN_%02d_N%02d_S1.csv" %(mags[erhan], np.abs(temps[erhan])) 
        data[erhan]['script1'] = pd.read_csv(script1)
        
        script2 = "dyn_data/THUN100_DYN_%02d_N%02d_S2.csv" %(mags[erhan], np.abs(temps[erhan]))
        data[erhan]['script2'] = pd.read_csv(script2)
        
        script3 = "dyn_data/THUN100_DYN_%02d_N%02d_S3.csv" %(mags[erhan], np.abs(temps[erhan]))  
        data[erhan]['script3'] = pd.read_csv(script3)
        
# Read modelocv.json file that contains parameters just for OCV data
# Such OCV0, OCVrel OCVeta Q etc.
with open('model_files/modelocv.json') as json_file:
    model = json.load(json_file)

# Create array for parameters g, RC1, M, R0dch, R0cha, and R1 
teta = np.arange(len(temps),dtype='float64')


# Set temperature values
model['Temps']=temps

# If it is firs time optimization then create spaces for parameter values in the model
if 'Q' not in model:
    not_first=False
    model['Q']=[0.0]*len(temps)
    model['eta']=[0.0]*len(temps)
    model['g']=[0.0]*len(temps)
    model['M']=[0.0]*len(temps)
    model['R0cha']=[0.0]*len(temps)
    model['R0dch']=[0.0]*len(temps)
    model['R1']=[0.0]*len(temps)
    model['RC1']=[0.0]*len(temps)
else:
    # Otherwise print exist parameters for information
    not_first=True
    print(model['RC1'])
    print(model['g'])
    print(model['M'])
    print(model['R0dch'])
    print(model['R0cha'])
    print(model['R1'])
    

        
# First we calculate actual capacity and coloumbic efficiency.
for erhan in range(len(temps)):
        
        s1disAh=data[erhan]['script1']['disAh'].iloc[-1]
        s2disAh=data[erhan]['script2']['disAh'].iloc[-1]
        s3disAh=data[erhan]['script3']['disAh'].iloc[-1]
        
        s1chgAh=data[erhan]['script1']['chgAh'].iloc[-1]
        s2chgAh=data[erhan]['script2']['chgAh'].iloc[-1]
        s3chgAh=data[erhan]['script3']['chgAh'].iloc[-1]
        
        
        # Obtain capacity from dynamic(UDDS) test data.
        model['Q'][erhan] =s1disAh + s2disAh - s1chgAh - s2chgAh
        
        # Obtain coulombic efficiency(eta) from dynamic(UDDS) test data.
        model['eta'][erhan] = (s1disAh + s2disAh+ s3disAh)/(s1chgAh + s2chgAh + s3chgAh)


# Boundary values for parameters. Lower bounds can't not be negative since
# they have physhical meaninngs and upper values are determined by experience
bounds=[[0.0]*len(teta),[1000, 1000, 1.0e-1, 1.0e-2, 1.0e-2, 1.0e-2]]

# Define callable function for least square estimation algorithm which
# takes parameter array as input and retur (simulation_values - true_values)
def fun(teta):
    fun.counter+=1

    # Copy paramters to simulate the model
    cell.RC1=teta[0]
    cell.g=teta[1]
    cell.M=teta[2]
    cell.R0dch=teta[3]
    cell.R0cha=teta[4]
    cell.R1=teta[5]
    cell.exp1 = np.exp(-cell.dt/cell.RC1)
    cell.one_exp1 = (1.0-cell.exp1)
    cell.dt_Q_3600 = cell.dt/(cell.Q * 3600.0) 
    cell.g_dt_Q = cell.g * cell.dt/cell.Q
    #print(temps[erhan])
    #cell.parameter_update(temps[erhan])
    # Calculate residual to return
    res = cell.fun([1.0,0,0],current)-voltage
    # Print some info after each iteration
    if  fun.counter % (len(teta)+1) == 0:
        rms = np.sqrt(np.mean(res**2))*1000
        print('Fun Evaluation Count=%d'%fun.counter)
        print('Root Mean Square Error=%fmV'%(rms))
        #print('Parameters')
        #print(teta)
    return res

# Options for least square algorithm for more details documentation 
# can be seen here
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
jac='cs'
method='trf'
machine_epsilon = np.finfo(float).eps
#xtol=ftol=gtol=machine_epsilon
xtol=ftol=gtol=1e-6
x_scale='jac'
loss='soft_l1'
f_scale=1
max_nfev=None
diff_step=None
tr_solver=None
tr_options={}
jac_sparsity=None
verbose=2
kwargs={}

# Define variables for optimisation results and rms error 
res={}
rmse = np.zeros(len(temps))


# Run estimation for all temperature values
for erhan in range(len(temps)):
    current = np.asarray(data[erhan]['script1']['current'])
    voltage = np.asarray(data[erhan]['script1']['voltage'])
    current = np.where( ( (current >= 0 ) & (current <= 0.1 ) ) , 0, current)
    current = np.where( ( (current <= 0 ) & (current >= -0.1 ) ) , 0, current)
    # Initial guess. It is important to choose reasonable valuess for fast covergence
    # If it is not first time estimation then use last estimation values as initial values
    if not_first:
        teta[0]= model['RC1'][erhan]
        teta[1]= model['g'][erhan]
        teta[2]= model['M'][erhan]
        teta[3]= model['R0dch'][erhan]
        teta[4]= model['R0cha'][erhan]
        teta[5]= model['R1'][erhan]
    else:
        # Otherwise use prepense initial values
        teta=[50, 20, 1.0e-2, 1.0e-3, 1.0e-3, 1.0e-3]
        model['RC1'][erhan]=1
        

    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0,use_OCVS=False)
    fun.counter=0

    # Start estimation for the specified temperature value
    print('Optimization Started for Temperature %02d °C' % temps[erhan])
    start2 = time.process_time()
    res[erhan] = least_squares(fun, teta, bounds=bounds,method=method, xtol=xtol, ftol=ftol, gtol=gtol, \
                        f_scale=f_scale, x_scale=x_scale,loss=loss, tr_solver=tr_solver, tr_options=tr_options, \
                        jac_sparsity=jac_sparsity, verbose=verbose, kwargs=kwargs)
    print(time.process_time() - start2)
    print('Optimization Finished!')
    print('Etimated Paramaters')
    print(res[erhan].x)
    
    # Save estimated parameters for related temperature value
    model['RC1'][erhan]=   res[erhan].x[0]
    model['g'][erhan]=     res[erhan].x[1]
    model['M'][erhan]=     res[erhan].x[2]
    model['R0dch'][erhan]= res[erhan].x[3]
    model['R0cha'][erhan]= res[erhan].x[4]
    model['R1'][erhan]=    res[erhan].x[5]


# Save the model file with estimated parameters
with open('DESC1Modeltest.json', 'w') as fp:
    json.dump(model, fp)
    
   
# Get stop time
now = datetime.now()
stop = now.strftime("%d/%m/%Y %H:%M:%S")

# Print Start and stop time
print("Start Time  =", start)	
print("Finish Time =", stop)	


# Plot estimated paramaters together
xfontsize = 12
yfontsize = 12
colors=['b', # blue
        'g', #green
        'r', #red
        'c', #cyan
        'm', #magenta
        'y',# yellow
        'orange', #orange
        'purple', #purple
        ]
temps=model['Temps']

for erhan in range(len(temps)):
    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0, use_OCVS=False)
    # Plot simulated output and cell data output
    current = np.asarray(data[erhan]['script1']['current'])
    voltage = np.asarray(data[erhan]['script1']['voltage'])
    time = np.asarray(data[erhan]['script1']['time'])/3600
    plt.figure()
    est_voltage = cell.fun([1.0,0,0],current)
    rmse[erhan]=np.sqrt(((est_voltage - voltage) ** 2).mean())
    plt.plot(time, est_voltage, color=colors[erhan], label='Simulation')
    plt.plot(time, voltage, label='True')
    plt.legend()
    plt.xlabel('Time(h)', fontsize = xfontsize, fontweight = 'bold')
    plt.ylabel('Voltage', fontsize = yfontsize, fontweight = 'bold')
    plt.title('True and Simulated Output Voltage for T=%02d °C, RMSE = %2.2fmV' % (temps[erhan],rmse[erhan]*1000))
    plt.show()

    # Save the plot
    plt.savefig('figures/estimations/param_est_temp_%02d.png' % temps[erhan], dpi=600, bbox_inches='tight')
    
    # Print rmse value
    print('Simulation with Estimated Parameters')
    print('RMS=%fmV'%(rmse[erhan]*1000)) 
    print('------------------------------------------------------------')


# Delete report file if exist
if os.path.exists("reports/estimations/param_est_report.txt"):
  os.remove("reports/estimations/param_est_report.txt")


# Plot rmse and save.
fig = plt.figure()
plt.plot(temps, rmse*1000, color=colors[0])
fig.suptitle('RMSE', fontsize = xfontsize, fontweight = 'bold')
plt.xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
plt.ylabel('RMSE of Output Voltage(mV)', fontsize = yfontsize, fontweight = 'bold')

# Tighten the plot and save
fig.tight_layout()
plt.savefig('figures/estimations/rmse_output_voltage_of_param_est.png', dpi=600, bbox_inches='tight')


# Create a report file and write results.
f = open("reports/estimations/param_est_report.txt", "x")
f.write('\n')
f.write('Parameters Estimation Results for DESC Model with 1-RC\n')
f.write('\n')

f.write('\n')
f.write('Temperature Values(°C):' + str(model['Temps']) + '\n')
f.write('\n')
f.write('RMSE Values(mV):' + str(rmse) + '\n')
f.write('\n')

fig,a =  plt.subplots(4,2)
a[0,0].plot(temps, model['M'], color=colors[0])
#a[0,0].set_title('Paramaeter M')
a[0,0].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[0,0].set_ylabel('M(mV)', fontsize = yfontsize, fontweight = 'bold')

a[0,1].plot(temps, model['g'], color=colors[1])
#a[0,1].set_title('Paramaeter \u03B3')
a[0,1].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[0,1].set_ylabel('\u03B3(unitless)', fontsize = yfontsize, fontweight = 'bold')

a[1,0].plot(temps, model['R0dch'], color=colors[2])
#a[0,1].set_title('Paramaeter R0dch')
a[1,0].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[1,0].set_ylabel('R0dch(m\u03A9)', fontsize = yfontsize, fontweight = 'bold')

a[1,1].plot(temps, model['R0cha'], color=colors[3])
#a[0,1].set_title('Paramaeter R0dch')
a[1,1].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[1,1].set_ylabel('R0dch(m\u03A9)', fontsize = yfontsize, fontweight = 'bold')

a[2,0].plot(temps, model['R1'], color=colors[4])
#a[0,1].set_title('Paramaeter R1')
a[2,0].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[2,0].set_ylabel('R1(m\u03A9)', fontsize = yfontsize, fontweight = 'bold')

a[2,1].plot(temps, model['RC1'], color=colors[5])
#a[0,1].set_title('Paramaeter RC1')
a[2,1].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[2,1].set_ylabel('RC1(s)', fontsize = yfontsize, fontweight = 'bold')

a[3,0].plot(temps,model['eta'], color=colors[6])
#a[0,1].set_title('Paramaeter eta')
a[3,0].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[3,0].set_ylabel('eta(unitless)', fontsize = yfontsize, fontweight = 'bold')

a[3,1].plot(temps,model['Q'], color=colors[7])
#a[0,1].set_title('Paramaeter \u03B3')
a[3,1].set_xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
a[3,1].set_ylabel('Q(Ah)', fontsize = yfontsize, fontweight = 'bold')

# Tighten the plot and save
fig.tight_layout()
plt.savefig('figures/estimations/estimated_parameters.png', dpi=600, bbox_inches='tight')


for erhan in range(len(temps)):
    f.write('True and Simulated Output Voltage RMSE=%fmV with Estimated Parameters at Temperature %02d°C'%(rmse[erhan]*1000, temps[erhan]) + '\n') 
f.write('\n')
f.write('Estimated Parameters\n')
f.write('Paramater M(mV):' + str(model['M']) + '\n')
#f.write('Paramater \u03B3(mV):' + str(model['g']) + '\n')
#f.write('Paramater R0dch(m\03a9):' + str(model['R0dch']) + '\n')
#f.write('Paramater R0cha(m\03a9):' + str(model['R0cha']) + '\n')
#.write('Paramater R1(m\03a9):' + str(model['R1']) + '\n')
f.write('Paramater g(mV):' + str(model['g']) + '\n')
f.write('Paramater R0dch(mohm):' + str(model['R0dch']) + '\n')
f.write('Paramater R0cha(mohm):' + str(model['R0cha']) + '\n')
f.write('Paramater R1(mohm):' + str(model['R1']) + '\n')
f.write('Paramater RC1(s):' + str(model['RC1']) + '\n')
f.write('Paramater eta(unitless):' + str(model['eta']) + '\n')
f.write('Paramater Cell Capacity Q(Ah):' + str(model['Q']) + '\n')
f.write('\n')
f.write("Start Time =" + start + '\n')
f.write("Stop Time =" +  stop + '\n')
f.write('\n')
f.close() # Close the file.


