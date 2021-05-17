# DESC_Model
 Double Enhanced Self Correcting Model


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

# For What is This Repo?
In this framework I have developed python code framework for my battery model and this framework also includes optimisation and parameter estimation scripts for my model.

https://data.mendeley.com/datasets/nbtrgh5frn/ here I have done similar thing in simulink for my research paper. Now moved the framework into Python.
Here https://github.com/batterysim/esctoolbox-python there is a similar framework for Enhanced Self Correcting Model(ESC). I used some part of codes here but the structure is totatly diferent.

DESC model is an electric equivalent circuit model for Lithium-Ion Batteries. This model can be modified for your own model and the same framework can be used. Also I have done lot of optimisation for my model (especialy with nonlinear least sqaures estimation) thus, one can use these optimisation scripts for different functions and it can be good start point.

# How to Start

## 1. OCV Curve Fitting
Run ocv.py file in the ocv_curve_fitting folder. This code provided from here https://github.com/batterysim/esctoolbox-python
ocv.py uses charge and discharge test data and takes its average to obtain rough OCV(Open Circuit Voltage) curves for all test temperatures and then does basic inter/extraploation for curve fitting.
That gives us good start point for further estimation.

## 2. Model Parameters Estimation from Dynmaic Data
Run DynamicProcess.py file. That start parameter estimation for DESC model with dynamic test data. Since I used 1-RC branch for my model there are 8 paramerters. These are Eta, Q, M, g, R1, RC1, R0cha, R0dch respectively. Eta and Q can be directly calculated from the test data and for the rest I run nonlinear least squares optimization from scipy.optimize.least_squares. There are many tricks and to lelarn about this library. Thus I recommend you to don't change settings at the beggining. For example tolerance values I made xtol=ftol=gtol=1e-6. This provides more aggreive optimisation but takes long time. After parameter estimation one can run ocv estimation and then again parameter estimation with lower tolerance values(like 1e-8). That provides better estimation results.

## 3. OCV Curves Estimation from Dynamic Data
At first step we have got OCV with basic averaging and inter/extrapolation. This gives as OCV curves for each temperature in a array with 201 elements.
So if we treat this arrays as parameter then we can do optimisation also on OCV curves. But this time we have 201 parameters instead of 6.
I used 2 methods for OCV estimation. First I have started with nonlinear least square(NLSQ) optimisation but I it generates some ripples on the curves these has no physcial meaning. Because as OCVs nature OCV should be decreased while discharging. Thus with NLSQ I have applied linear inter/extrapolation to smooth this curves. It gives reasonable results somehow.

Then I have discovered sequantional Least Squares Programming(SLSQP) library in python which was what I needed. Since I am not so familiar with optimisation I didn't know constrained optimisation. SLSQP provides you constrained optimisation so I used constrain as OCV[n]>=OCV[n-1]. That means OCV values cannot be bigger than previous values. SLSQP doesn't require interpolation. Thus, it is more convinient and gives better results then nice OCV curves then NLSQ. How ever I kept both files for comparision.

## Steps 2-3 Can be Done Sequentially
After first you have estimated OCV curves and parameters then one can run step 2-3 sequentially with lower tolerance values. That provide more aggresive estimation but takes longer time. For example for OCV estimation once I had to wait more than 10 hours.


## Simulation and Plot Data
DESC_Simulation.py can be used to simulate model with estimated parameters.
Plot_Dynamic_Data.py can be used to plot dynamic test data
Plot_OCV_Data.py can be used to plot OCV test data
