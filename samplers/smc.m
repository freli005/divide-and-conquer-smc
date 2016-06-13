function [X,lZ,logW,ess_log] = smc(X,Jv0,Jh0,Jv1,Jh1,alpha_vec,par)
% Input:
%   X, M x M x N: initial samples
%   Jv0, Jh0, M x M (or scalar): initial edge interactions
%   Jv1, Jh1, M x M: final edge interactions
%   alpha_vec, 1 x (n+1): Annealing sequence, where n is the number of annealing
%                        steps; alpha_vec(1) = 0, alpha_vec(n+1) = 1
%   par.resampling: 0 = off (AIS), 1 = multinomial, 2 = stratified, 3 = systematic
%      .ESS_threshold: Value in [0,1], determining when to resample
%
% Output:
%   X, M x M x N: final samples
%   lZ, scalar: estimate of log[Z(1)/Z(0)]
%   logW, N x 1: final unnormalised log-weights
%   ess_log, 1 x (n+1): ESS values at each iteration

if(alpha_vec(1) ~= 0)
    error(message('smc:alphaInitNotZero'));
end

% Analytics
ess_log = zeros(1,length(alpha_vec));
ess_log(1) = 1;

% Setup
[nx,ny,N] = size(X);
logW = zeros(N,1);
lZ = 0;
Jvt = Jv1-Jv0; % Compute incremental edge potentials / J = J0 + alpha*Jt
Jht = Jh1-Jh0;

numAlpha = length(alpha_vec)-1;
for iAlpha = 1:numAlpha
    %%% Propagate from alpha --> alpha_new
    if(iAlpha > 1) % Don't propagate on the first iteration
        X = gibbsgrid(X, alpha_vec(iAlpha)*Jvt + Jv0, alpha_vec(iAlpha)*Jht + Jh0);
    end
    
    %%% Weight computation    
    % Computes sum_{ij} J_{ij}*x_i*x_j for every page in X
    S = sum_diff_x(X, Jvt, Jht);    
    % Compute new weights
    logW =  logW + (alpha_vec(iAlpha+1)-alpha_vec(iAlpha))*S;
    
    %%% Compute ESS / resample if needed
    maxlW = max( logW );
    w = exp( logW - maxlW );
    W = w/sum(w);
    ess = 1/(N*sum(W.^2));
    % Store for analysis
    ess_log(iAlpha+1) = ess;
    
    %%% Resample if needed
    if(iAlpha < numAlpha && par.resampling > 0 && ess < par.ESS_threshold) % Never resample at last iterations
        ind = resampling(W,par.resampling);
        X = X(:,:,ind);
        % Update normalizing constant at this iteration
        lZ = lZ + maxlW + log( sum( w ) ) - log( N );
        % Reset weights
        logW = zeros(N,1);
    end
end