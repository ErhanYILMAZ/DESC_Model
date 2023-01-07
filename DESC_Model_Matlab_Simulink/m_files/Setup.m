% This m file contains basic code scripts to provide UDDS data from m
% files to simulink models for parameter estimation. As wel as to get rmse
% and mae values.

% Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
% This work is licensed under CC BY-SA (Attribution-ShareAlike)
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.

%% Calculate RMS and MA Error Values for Identified Parameters

rms_errorr= sqrt(mean(diff.data(1,:).^2));
ma_errorr= mean(abs(diff.data(1,:)));


%% Load UDDS Test Data and Create Time Series.
Current=THUN_UDDS_Data(1).script1.current;
Voltage=THUN_UDDS_Data(1).script1.voltage;
t=(0:length(Current)-1)';
%%
