#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
This file is licensed under CC BY-SA (Attribution-ShareAlike)
It is provided "as is", without express or implied warranty, for educational and informational purposes only.

Python file simulate Double Enhanced Self Correcting(DESC) model with 1-RC branch

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
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from DESC_1RC import DESC_1RC
from datetime import datetime

now = datetime.now()
start = now.strftime("%d/%m/%Y %H:%M:%S")


temps=[-5, 5, 15, 25, 35, 45]
mags= [20, 20, 40, 50, 50, 50]


# Read modelocv.json file that contains parameters just for OCV data
# Such OCV0, OCVrel OCVeta Q etc.
with open('model_files/DESC1Model_SLSQP.json') as json_file:
    model = json.load(json_file)

# Define variable for rmse value
rmse = np.zeros(len(temps))

# Define font size for labels
xfontsize = 12
yfontsize = 12

for erhan in range(len(temps)):

    # Read UDDS dynamic data for specified temperature value
    if(temps[erhan]>0):
        script1  = "dyn_data/THUN100_DYN_%02d_P%02d_S1.csv" %(mags[erhan], temps[erhan]) 
        data = pd.read_csv(script1)
    else:
        script1  = "dyn_data/THUN100_DYN_%02d_N%02d_S1.csv" %(mags[erhan], np.abs(temps[erhan])) 
        data= pd.read_csv(script1)
    
    # Get voltage and current values.
    current = np.asarray(data['current'])
    voltage = np.asarray(data['voltage'])
    time = np.asarray(data['time'])/3600
    # Create model instance and initiliaze for the specified temperature.
    cell = DESC_1RC(model=model , temp=temps[erhan], dt=1.0,use_OCVS=True)

    # Plot simulated output and cell data output
    plt.figure()
    est_voltage = cell.fun([1.0,0.0,0.0],current)
    rmse[erhan]=np.sqrt(((voltage - est_voltage) ** 2).mean())
    plt.plot(time, est_voltage)
    plt.plot(time, voltage)
    plt.xlabel('Time(h)', fontsize = xfontsize, fontweight = 'bold')
    plt.ylabel('Voltage', fontsize = xfontsize, fontweight = 'bold')
    plt.title('Estimated Output Voltage for T=%02d, RMSE = %2.2fmV' % (temps[erhan],rmse[erhan]*1000))
    plt.show()
    plt.savefig('figures/simulations/Simulation_T_%02d_RMSE_%2.2fmV.png' % (temps[erhan],rmse[erhan]*1000), dpi=600, bbox_inches='tight')
    # Print rms error.
    print('RMS=%fmV'%(rmse[erhan]*1000))
    print('------------------------------------------------------------')

# Get stop time
now = datetime.now()
stop = now.strftime("%d/%m/%Y %H:%M:%S")

# Delete report file if exist
if os.path.exists("reports/simulations/simulation_report.txt"):
  os.remove("reports/simulations/simulation_report.txt")

# Create a report file and write results.
f = open("reports/simulations/simulation_report.txt", "a")
f.write('\r\n')
f.write('Simulation Results for DESC Model with 1-RC\n')
f.write('\r\n')
for erhan in range(len(temps)):
    f.write('Simulated Output Voltage RMSE=%2.2fmV at %02d°C Temp\n' % (rmse[erhan]*1000, temps[erhan])) 
f.write('\r\n')
f.write("Start Time =" + start)
f.write("\nStop Time =" +  stop)
f.write('\n')
f.close() # Close the file.

