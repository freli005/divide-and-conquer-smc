% Runs a standard SMC sampler on the Ising model
clear;
addpath('./utils');
addpath('./samplers');
[Jv,Jh] = setup_model(6); % M = 2^6 = 64

%% Set some parameters for the D&C-SMC sampler
N = 64;
% method = 'dcsir';
% method = 'dcsmc (ann)';
% method = 'dcsmc (mix)';
% method = 'dcsmc (ann+mix)';
method = 'dcsmc (ann+lw)';

switch method
    case 'dcsir'
        par.merge_method = 0;       % Gives SIR style resampling (merge step)
        par.CESS_threshold = 0;     % Results in no annealing, the adaptive procedure will always select \alpha = [0,1]        
        par.ESS_threshold = 0.5;    % Not used by D&C-SIR, really, but the dcsmc function looks at this value to compute some statistics so it has to be set to something!
    case 'dcsmc (ann)'
        par.merge_method = 0;       % Gives SIR style resampling (merge step)
        % The three settings below affect the "SMC sampler" (i.e. tempering) phases of the algorithm
        par.resampling = 3;         % 0 = off (AIS), 1 = multinomial, 2 = systematic, 3 = stratified
        par.ESS_threshold = 0.5;    % Value in [0,1], determining when to resample
        par.CESS_threshold = 0.995; % Value in [0,1], CESS target value used when adapting alpha
    case 'dcsmc (mix)'
        par.merge_method = 2;       % Mixture sampling
        par.ESS_merge_threshold = 0; % Always merge at alpha* = 1
        par.merge_subsampling = 2;  % Determines how the N^2 system is subsampled (2 = stratified on marginal/conditional, see smc_merge.m for more details)
        par.CESS_threshold = 0;     % Results in no annealing, the adaptive procedure will always select \alpha = [0,1]        
        par.ESS_threshold = 0.5;    % Not used by D&C-SMC (mix), really, but the dcsmc function looks at this value to compute some statistics so it has to be set to something!
    case 'dcsmc (ann+mix)'
        par.merge_method = 2;       % Mixture sampling (adaptation based on marginal CESS)
        par.ESS_merge_threshold = 0.95; % (C)ESS thresholf for adapting alpha*
        par.merge_subsampling = 2;  % Determines how the N^2 system is subsampled (2 = stratified on marginal/conditional, see smc_merge.m for more details)
        % The three settings below affect the "SMC sampler" (i.e. tempering) phases of the algorithm
        par.resampling = 3;         % 0 = off (AIS), 1 = multinomial, 2 = systematic, 3 = stratified
        par.ESS_threshold = 0.5;    % Value in [0,1], determining when to resample
        par.CESS_threshold = 0.995; % Value in [0,1], CESS target value used when adapting alpha
    case 'dcsmc (ann+lw)'
        par.merge_method = 3;       % Lightweight mixture sampling (adaptation based on ESS of the mN system)        
        par.lw_m = 32;              % subsamples the N^2 population "blindly" to m*N samples before adapting alpha*
        par.ESS_merge_threshold = 16; % Select alpha* so that the ESS of the mN particle system is  par.ESS_merge_threshold*N (note that we require par.ESS_merge_threshold <= par.lw_m)        
        % The three settings below affect the "SMC sampler" (i.e. tempering) phases of the algorithm
        par.resampling = 3;         % 0 = off (AIS), 1 = multinomial, 2 = systematic, 3 = stratified
        par.ESS_threshold = 0.5;    % Value in [0,1], determining when to resample
        par.CESS_threshold = 0.995; % Value in [0,1], CESS target value used when adapting alpha                
        % Do not change these for case (lw) - will result in error message
        par.merge_subsampling = 1;  % Determines how the m*N particle system i subsampled. Has to be set to 1 for case (lw) ==> direct stratified subsampling
        par.merge_exponent = 1;     % Determines a scaling for ESS_merge_theshold. Has to be set to 1 for case (lw) ==> we measure ESS of m*N particle systems
    otherwise
        error('Unknown method');
end

par.printOn = 0; % Turn on to show progress of algorithm


%%% Some default shared settings below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.CESS_tolerance = 1e-5; % Tolerance for the newton/bisection method used to adapt alpha

% We experimented with using a smoothed version of the reintroduction of
% the missing edges between subpopulations. That is, instead of raising the
% annealing parameter of all edges within a subpopulation to 1 before
% proceeding to the next merge step, we kept the parameter at a value < 1
% for any MRF edges close to the boundary of each subpopulation lattice.
% However, this gave rise to a tuning parameter which was difficult to
% reason about and, thus, we did not use this idea for the results
% presented in the paper (par.mask_order = 1). Setting par.mask_order > 1
% (integer) results in smoothing of the edges.
par.mask_order = 1; % 1 = no smoothing, as used in paper


%% Run D&C-SMC sampler
tic;
[X, W, lZ, totMCMC, totRes, alphaStar] = dcsmc(Jv, Jh, N, [], par);
cpu_time = toc;

% Expected energy
engy = -sum_diff_x(X,Jv,Jh); % Energy
Ehat = (engy')*W(:);

%% Print results
fprintf('\n\n\nD&C-SMC Sampler finised in %2.2f seconds.\n',cpu_time);
fprintf('On average, %2.2f MCMC updates and %2.2f resampling steps were carried out for each site.\n',mean(mean(sum(totMCMC,3))),mean(mean(sum(totRes,3))));
fprintf('The estimate of the log-normalising-constant is: %2.2f\n',lZ);
fprintf('The estimate of the expected energy is: %2.2f\n',Ehat);