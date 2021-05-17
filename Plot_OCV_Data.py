#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
This file is licensed under CC BY-SA (Attribution-ShareAlike)
It is provided "as is", without express or implied warranty, for educational and informational purposes only.

This python script plots OCV data collected for
ThunderSky-Winston-LIFEPO4-100Ah-WIDE cell according to test descriptions below in the files.

http://mocha-java.uccs.edu/BMS1/index.html
http://mocha-java.uccs.edu/ECE5710/index.html
http://mocha-java.uccs.edu/ECE5710/ECE5710-Notes02.pdf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

temps=[-5, 5, 15, 25, 35, 45]


for erhan in range(len(temps)):

    for erhan2 in range(4):
        if(temps[erhan]>0):
            script  = "ocv_data/THUN100_OCV_P%02d_S%1d.csv" % (temps[erhan], erhan2+1)
        else:
            script  = "ocv_data/THUN100_OCV_N%02d_S%1d.csv" % (np.abs(temps[erhan]), erhan2+1)
            
        data = pd.read_csv(script)
        current = -np.asarray(data['Current'])
        voltage = np.asarray(data['VCell1'])
        time    =  np.arange(start=0, stop=voltage.size, step=1)/360.0

        fig,a =  plt.subplots(2,1)
        x = np.arange(1,5)
        a[0].plot(time, current)
        a[0].set_title('OCV Script %1d T=%02d' % (erhan2+1, temps[erhan]))
        #a[0].set_xlabel('Time(h)')
        a[0].set_ylabel('Current(A)')
        a[1].plot(time, voltage)
        #a[1].set_title('OCV Script %1d T=%02d' % (erhan2+1, temps[erhan]))
        a[1].set_xlabel('Time(h)')
        a[1].set_ylabel('Voltage(V)')
        plt.show()


