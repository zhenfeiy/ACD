% It will first simulate a ACD model and then fit it using the fitting function

clear;

addpath('Func');

nr=1000;    % No. of observations in the simulation

dist='exp';  % dist='exp' or 'weibull';

Coeff.w=.1893; % constant in expected duration (psi)
Coeff.q=.1422; % Coeff at duration in t-1 (alpha)
Coeff.p=.6691; % Coeff at expected duration in t-1 (beta)
Coeff.y=.8; % just for weibull dist

q=size(Coeff.q,2);
p=size(Coeff.p,2);

simulDur=ACD_Simul(nr,Coeff,q,p,dist);  % Simulation

rmpath('Func');
