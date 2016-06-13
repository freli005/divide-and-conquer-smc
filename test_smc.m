% Runs a standard SMC sampler on the Ising model
clear;
addpath('./utils');
addpath('./samplers');
[Jv,Jh] = setup_model(6); % M = 2^6 = 64

%% Set some parameters for the SMC sampler
N = 64;
par.resampling = 3; % 0 = off (AIS), 1 = multinomial, 2 = systematic, 3 = stratified
par.ESS_threshold = 0.5; % Value in [0,1], determining when to resample
par.CESS_threshold = 0.995; % Value in [0,1], CESS target value used when adapting alpha
par.CESS_tolerance = 1e-5; % Tolerance for the newton/bisection method used to adapt alpha

par.printOn = 0; % Turn on to show progress of algorithm

%% Run SMC sampler
tic;
[X, W, lZ, alpha_log, ess_log] = smcsampler(Jv, Jh, N, [], par);
cpu_time = toc;
% The vector beta_log contains all annealing temperatures, including [0,1].
% Hence, the total number of applications of the MCMC kernel is:
totMCMC = length(alpha_log)-2;
% Count the number of resampling steps that are taken
totRes = sum(ess_log(1:end-1) <= par.ESS_threshold);
% Expected energy
engy = -sum_diff_x(X,Jv,Jh); % Energy
Ehat = (engy')*W(:);

%% Print results
clc
fprintf('SMC Sampler finised in %2.2f seconds (%i iterations).\n',cpu_time,totMCMC);
fprintf('In total %i resampling steps were carried out.\n',totRes);
fprintf('The estimate of the log-normalising-constant is: %2.2f\n',lZ);
fprintf('The estimate of the expected energy is: %2.2f\n',Ehat);

figure(1);
plot(alpha_log);
xlabel('Iteration');
ylabel('\alpha');