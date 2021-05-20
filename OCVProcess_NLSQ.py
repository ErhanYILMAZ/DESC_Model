#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
This file is licensed under CC BY-SA (Attribution-ShareAlike)
It is provided "as is", without express or implied warranty, for educational and informational purposes only.

Python file for Double Enhanced Self Correcting(DESC) model with 1-RC branch parameter estimation.
This file does open circuit voltage (OCV) curves estimation for model from dynamic test data and uses 
nonlinear least squares(NLSQ) optimisation for parameter estimation.

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
with open('model_files/DESC1Model.json','r') as json_file:
    model = json.load(json_file)


# Set temperature values
model['Temps']=temps

# Define callable function for least square estimation algorithm which
# takes parameter array as input and retur (simulation_values - true_values)
def fun(teta):
    fun.counter+=1
    # Copy paramters to simulate the model
    cell.ocv=teta
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
fun.counter=0


# Options for least square algorithm for more details documentation 
# can be seen here
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
# https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#least-squares-minimization-least-squares

jac='cs'
method='trf'
xtol=ftol=gtol=1e-2
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

# Define variables for optimisation results and rms error 
OCVs_Raw={}
res={}
rmse = np.zeros(len(temps))

for erhan in range(len(temps)):
    # Get current and voltage values for the specified test temperature
    current = np.asarray(data[erhan]['script1']['current'])
    voltage = np.asarray(data[erhan]['script1']['voltage'])
    current = np.where( ( (current >= 0 ) & (current <= 0.1 ) ) , 0, current)
    current = np.where( ( (current <= 0 ) & (current >= -0.1 ) ) , 0, current)
    
    # We initialize our Double Enhanced Self Correcting Model wit 1-RC branch
    # at specified test emperature and sample period is 1s

    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0, use_OCVS=False)
    fun.counter=0

    teta = cell.ocv
        
    # Boundary values for parameters. lower bounds can't not be negative since
    # they have physhical meaninngs and upper values are determined by experience
    bounds=[[2.6]*len(teta),[3.8]*len(teta)]

    # Start estimation for the specified temperature value
    print('Optimization Started for Temperature %02d °C' % temps[erhan])
    start2 = time.process_time()
    res[erhan] = least_squares(fun, teta, bounds=bounds,method=method, xtol=xtol, ftol=ftol, gtol=gtol, \
                        f_scale=f_scale, x_scale=x_scale,loss=loss, tr_solver=tr_solver, tr_options=tr_options, \
                         diff_step=diff_step,max_nfev=max_nfev, jac_sparsity=jac_sparsity,  verbose=verbose, kwargs=kwargs)
    print(time.process_time() - start2)
    print('Optimization Finished!')
    print('Etimated Paramaters')
    print(res[erhan].x)
    
    # Save estimated parameters for related temperature value
    OCVs_Raw[erhan] = res[erhan].x

    
# NLS Optimisation creates ripples on the OCV curves those don't effect models behavior
# but has no physical meaning. Thus we have used polynomial fitting here to get smooth ocv curves.
# Thus smoothing might change rms error sligthy but this is trivial.

soc=np.linspace(0,1,len(model['SOC']))
OCVs={}
for erhan in range(len(temps)):
    mymodel = np.poly1d(np.polyfit(soc, OCVs_Raw[erhan], 500))
    OCVs[erhan] = mymodel(soc)
    plt.figure()
    plt.plot(soc*100, OCVs_Raw[erhan], color=colors[erhan], label='NLSQ + Interp')
    plt.plot(soc*100, OCVs[erhan], label='NLSQ')
    plt.legend()
    plt.xlabel('SoC(%)', fontsize = xfontsize, fontweight = 'bold')
    plt.ylabel('Voltage(V)', fontsize = yfontsize, fontweight = 'bold')
    plt.title('Estimated SOC Curve at T=%02d °C' % (temps[erhan]))
    plt.show()
    
    #Save the plot
    plt.savefig('figures/OCV_NLSQP_est/estimated_OCV_with_NLSQ_at_temp_%02d.png' % temps[erhan], dpi=600, bbox_inches='tight')
  
# get stop time an format in as in proper string
now = datetime.now()
stop = now.strftime("%d/%m/%Y %H:%M:%S")

# Print start and stop time
print("Start Time  =", start)	
print("Finish Time =", stop)	

# After OCV estimation for each temperature, OCV0 and OCVrel can be
# obtained below by the least squares estimation. This gives moreles same
# accuracy as using OCVs curves for each temperature. Also since it is
# only OCV0 and OCVrel it stores less space in the model files.

Y = np.zeros([len(temps), len(model['SOC'])])   # initialize rawocv array

for erhan in range(len(temps)):
    Y[erhan] = OCVs[erhan]

A = np.ones([len(temps), 2])
A[:, 1] = temps
X = np.dot(np.linalg.pinv(A) , Y)

# Store estimated and calculated ocv curves in the model
model['OCV0'] = list(X[0])
model['OCVrel'] = list(X[1])

model['OCV_N05'] = list(OCVs[0])
model['OCV_P05'] = list(OCVs[1])
model['OCV_P15'] = list(OCVs[2])
model['OCV_P25'] = list(OCVs[3])
model['OCV_P35'] = list(OCVs[4])
model['OCV_P45'] = list(OCVs[5])

# Save the model file with name stamp NLSQ that specifies
# non linear least squares method used for estimation
with open('model_files/DESC1Model_NLSQ.json', 'w') as fp:
    json.dump(model, fp)


# Get stop time
now = datetime.now()
stop = now.strftime("%d/%m/%Y %H:%M:%S")

# Print Start and stop time
print("Start Time  =", start)	
print("Finish Time =", stop)	


# Delete report file if exist
if os.path.exists("reports/estimations/OCV_NLSQ_est_report.txt"):
  os.remove("reports/estimations/OCV_NLSQ_est_report.txt")


# Create a report file and write results.
f = open("reports/estimations/OCV_NLSQ_est_report.txt", "x")
f.write('\n')
f.write('OCV Curves Fitting with (NLSQ) Results  for DESC Model with 1-RC\n')
f.write('\n')

for erhan in range(len(temps)):
    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0, use_OCVS=True)
    # Get current, voltage, and time values for the specified test temperature
    current = np.asarray(data[erhan]['script1']['current'])
    voltage = np.asarray(data[erhan]['script1']['voltage'])
    time = np.asarray(data[erhan]['script1']['time'])/3600
    current = np.where( ( (current >= 0 ) & (current <= 0.1 ) ) , 0, current)
    current = np.where( ( (current <= 0 ) & (current >= -0.1 ) ) , 0, current)
    
    # Plot simulated output and cell data output
    plt.figure()
    est_voltage = cell.fun([1.0,0,0],current)
    rmse[erhan]=np.sqrt(((est_voltage - voltage) ** 2).mean())
    plt.plot(time, est_voltage, color=colors[erhan], label='Simulation')
    plt.plot(time, voltage, label='True')
    plt.legend()
    plt.xlabel('Time(h)', fontsize = xfontsize, fontweight = 'bold')
    plt.ylabel('Voltage', fontsize = yfontsize, fontweight = 'bold')
    plt.title('True and Simulated Output with Estimated OCV T=%02d °C, RMSE = %2.2fmV' % (temps[erhan],rmse[erhan]*1000))
    plt.show()
    
    # Save the plot
    plt.savefig('figures/OCV_NLSQP_est/OCV_NLSQ_est_with_OCVs_temp_%02d.png' % temps[erhan], dpi=600, bbox_inches='tight')
    
    # Print rmse value
    print('Simulation with OCV Curves')
    print('RMS=%fmV'%(rmse[erhan]*1000))
    print('------------------------------------------------------------')


# Plot rmse and save.
fig = plt.figure()
plt.plot(temps, rmse*1000, color=colors[0])
fig.suptitle('RMSE', fontsize = xfontsize, fontweight = 'bold')
plt.xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
plt.ylabel('RMSE of Output Voltage(mV)', fontsize = yfontsize, fontweight = 'bold')

# Tighten the plot and save
fig.tight_layout()
plt.savefig('figures/OCV_NLSQP_est/rmse_output_voltage_of_NLSQ_OCV_est.png', dpi=600, bbox_inches='tight')


f.write('\n')
for erhan in range(len(temps)):
    f.write('True and SImulated Output Voltage RMSE=%fmV with Estimated OCV Curve at Temperature %02d°C'%(rmse[erhan]*1000, temps[erhan]) + '\n')
f.write('\n')
f.write('Temperature Values(°C):' + str(model['Temps']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature -05°C:' + str(model['OCV_N05']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature +05°C:' + str(model['OCV_P05']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature +15°C:' + str(model['OCV_P15']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature +25°C:' + str(model['OCV_P25']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature +35°C:' + str(model['OCV_P35']) + '\n')
f.write('\n')
f.write('OCV Values at Temperature +45°C:' + str(model['OCV_P45']) + '\n')
f.write('\n')

for erhan in range(len(temps)):
    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0, use_OCVS=False)
    # Get current, voltage, and time values for the specified test temperature
    current = np.asarray(data[erhan]['script1']['current'])
    voltage = np.asarray(data[erhan]['script1']['voltage'])
    time = np.asarray(data[erhan]['script1']['time'])/3600
    current = np.where( ( (current >= 0 ) & (current <= 0.1 ) ) , 0, current)
    current = np.where( ( (current <= 0 ) & (current >= -0.1 ) ) , 0, current)
    # Plot simulated output and cell data output
    plt.figure()
    est_voltage = cell.fun([1.0,0,0],current)
    rmse[erhan]=np.sqrt(((est_voltage - voltage) ** 2).mean())
    plt.plot(time, est_voltage, color=colors[erhan], label='Simulation')
    plt.plot(time, voltage, label='True')
    plt.legend()
    plt.xlabel('Time(h)', fontsize = xfontsize, fontweight = 'bold')
    plt.ylabel('Voltage', fontsize = yfontsize, fontweight = 'bold')
    plt.title('True and Simulated Output with Estimated OCV0 and OCVrel T=%02d °C, RMSE = %2.2fmV' % (temps[erhan],rmse[erhan]*1000))
    plt.show()

    # Save the plot
    plt.savefig('figures/OCV_NLSQP_est/OCV_NLSQ_est_with_OCV0_and_OCVrel_temp_%02d.png' % temps[erhan], dpi=600, bbox_inches='tight')
    
    # Print rmse value
    print('Simulation with OCV0 and OCVrel')
    print('RMS=%fmV'%(rmse[erhan]*1000))
    print('------------------------------------------------------------')

# Plot rmse and save.
fig = plt.figure()
plt.plot(temps, rmse*1000, color=colors[0])
fig.suptitle('RMSE', fontsize = xfontsize, fontweight = 'bold')
plt.xlabel('Temp(°C)', fontsize = xfontsize, fontweight = 'bold')
plt.ylabel('RMSE of Output Voltage(mV)', fontsize = yfontsize, fontweight = 'bold')

# Tighten the plot and save
fig.tight_layout()
plt.savefig('figures/OCV_NLSQP_est/rmse_output_voltage_of_NLSQ_OCV0_OCVrel_est.png', dpi=600, bbox_inches='tight')



f.write('\n')
for erhan in range(len(temps)):
    f.write('True and Simulated Output Voltage RMSE=%fmV with Estimated OCV0 and OCVrel at Temperature %02d°C'%(rmse[erhan]*1000, temps[erhan]) + '\n')
f.write('\n')
f.write('Temperature Values(°C):' + str(model['Temps']) + '\n')
f.write('\n')
f.write('OCV0 Values:' + str(model['OCV0']) + '\n')
f.write('\n')
f.write('OCVrel Values:' + str(model['OCVrel']) + '\n')
f.write('\n')
f.write("Start Time =" + start + '\n')
f.write("Stop Time =" +  stop + '\n')
f.write('\n')
f.close() # Close the file.



