% This will fisrt simulate a ACD (1,1)
% and then plot the likelihood of the model
clear;

addpath('Func');

nr=1000;            
n_simulations=5000; % How many simulations at likelihood plot

Coeff.q=.2; % Coeff at duration in t-1 (alpha)
Coeff.p=.7; % Coeff at expected duration in t-1(beta)

Coeff.w=.1;

simulDur=ACD_Simul(nr,Coeff,1,1,'exp');  % Simulation (exp or weibull)

plotLikelihoodFunct(simulDur,n_simulations); % This is the function for plotting the log likelihoods

rmpath('Func');
