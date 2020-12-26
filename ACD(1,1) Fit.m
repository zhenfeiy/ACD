% It will first load some duration data and then fit a ACD(1,1) by using the fitting function

clear;
load addur.mat;
addpath('Func');
dur = addur;

dist='exp';   % dist='exp' or 'weibull';

stdMethod=1;   % standard error calculation

q=1;
p=1;
x=dur;

[specOut]=ACD_Fit(x,dist,q,p,stdMethod);    % Fitting

plot([specOut.h]);
title('Duration Simulation and Modelling');
legend('Fitted Duration', 'Real Duration');
xlabel('Observations');
ylabel('Durations');

rmpath('Func');
