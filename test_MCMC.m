% Runs a single flip MH or blocked Gibbs sampler on the Ising model
clear;
addpath('./utils');
addpath('./samplers');
[Jv,Jh] = setup_model(6); % M = 2^6 = 64

%% Set some parameters for the MCMC sampler
numMCMC = 2^10;
par.blocksize = 1; % 1 results in single flip MH and > 1 results in blocked Gibbs
par.plotOn = 0;  % Plot state at each iteration; for B > 1, set plotOn = 2 to plot intermediate block updates
par.printOn = 0; % Turn on to show progress of algorithm

%% Run MCMC sampler
tic;
[X, Ehat] = gibbssampler(Jv, Jh, numMCMC, 1, par);
cpu_time = toc;

%% Print results
clc
fprintf('MCMC Sampler finised in %2.2f seconds (%i iterations).\n',cpu_time,numMCMC);
fprintf('The estimate of the expected energy is: %2.2f\n',Ehat);