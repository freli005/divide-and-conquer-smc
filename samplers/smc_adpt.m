function [X,lZ,logW,alpha_log,ess_log] = smc_adpt(X,Jv0,Jh0,Jv1,Jh1,alpha_init,par)
% Input:
%   X, M x M x N: initial samples
%   Jv0, Jh0, M x M (or scalar): initial edge interactions
%   Jv1, Jh1, M x M: final edge interactions
%   par.resampling: 0 = off (AIS), 1 = multinomial, 2 = stratified, 3 = systematic
%      .ESS_threshold: Value in [0,1], determining when to resample
%      .CESS_threshold: Value in [0,1], CESS target value used when adapting alpha
%
% Output:
%   X, M x M x N: final samples
%   lZ, scalar: estimate of log[Z(1)/Z(0)]
%   logW, N x 1: log-importance-weights (unnormalised)
%   alpha_log, 1 x (n+1): increments in the annealing schedule
%   ess_log, 1 x (n+1): ESS values at each iteration


% Setup
alpha_log = zeros(1,1000); % Store the alpha-values that we generate for analysis purposes
alpha = alpha_init; % Current value of alpha
alpha_log(1) = alpha;

ess_log = zeros(1,1000);
ess_log(1) = 1;
k = 1;

[nx,ny,N] = size(X);
logW = zeros(N,1);
lZ = 0;
Jvt = Jv1-Jv0; % Compute incremental edge potentials / J = J0 + alpha*Jt
Jht = Jh1-Jh0;

alpha_increment = 0.01; % Used to set alpha_init at first iteration of the adaptation
this_is_smc = strcmpi(par.sampler,'smc');
while(alpha < 1)
    if(par.printOn && this_is_smc)
        fprintf('alpha = %2.6f\n',alpha);
    end
    
    %%% Propagate from alpha --> alpha_new
    if(k > 1) % Don't propagate on the first iteration
         X = gibbsgrid(X, alpha*Jvt + Jv0, alpha*Jht + Jh0);
    end
    
    k = k+1;
    
    %%% Step size selection / weight computation
    % Computes sum_{ij} J_{ij}*x_i*x_j for every page in X
    S = sum_diff_x(X, Jvt, Jht);    
    % Adapt the annealing schedule based on CESS    
    alpha_init = min( (1-alpha)/2, alpha_increment );
    [logW, alpha_increment, ess] = adpt_alpha(S, logW, alpha_init, 1-alpha, par.CESS_threshold, par.CESS_tolerance); % Adapt alpha and compute new weights based on CESS criteria
    alpha = alpha + alpha_increment;    
    % Store for analysis
    alpha_log(k) = alpha_increment;
    ess_log(k) = ess;
    
    %%% Resample if needed
    if(alpha < 1 && par.resampling > 0 && ess < par.ESS_threshold) % Never resample at last iterations
        maxlW = max( logW );
        w = exp( logW - maxlW );
        W = w/sum(w);

        ind = resampling(W,par.resampling);
        X = X(:,:,ind);
        % Update normalizing constant at this iteration
        lZ = lZ + maxlW + log( sum( w ) ) - log( N );
        % Reset weights
        logW = zeros(N,1);
    end
end

% Truncate
alpha_log = alpha_log(1:k);
ess_log = ess_log(1:k);