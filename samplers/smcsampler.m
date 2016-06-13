function [X,W,lZ,alpha_log,ess_log] = smcsampler(Jv, Jh, N, alpha, par)
% Only works for square models, MxM. Can handle periodic or free boundary
% conditions.
%
% Jv, M x M: vertical interactions, Jv(i,j) is the interaction from X(i,j)
%            to X(i+1,j). The M'th row of Jv containts the boundary
%            interactions between X(M,j) and X(1,j). If omitted, this is
%            assumed to be zero (free boundary condition).
% Jh, M x M: horizontal interactions, Jh(i,j) is the interaction from
%            X(i,j) to X(i,j+1). The M'th column of Jh contains the
%            boundary interactions. If omitted, this is assumed to be zero
%            (free boundary condition).
% N, integer: number of particles
% alpha, 1 x (n+1): Annealing sequence, where n is the number of annealing
%                  steps. alpha(1) = 0, alpha(n+1) = 1
%       []: Adaptive annealing
% par.resampling: 0 = off (AIS), 1 = multinomial, 2 = stratified, 3 = systematic
%    .ESS_threshold: Value in [0,1], determining when to resample
%    .CESS_threshold: Value in [0,1], CESS target value used when adapting alpha

par.sampler = 'smc';
adaptive = numel(alpha)==0;
M = size(Jv,2); % Size of grid

% For free boundary conditions, assume periodic with zero interaction
if(size(Jv,1) < M)
    Jv = [Jv ; zeros(1,M)];
    Jh = [Jh zeros(M,1)];
end

nx = M;
ny = M;


% "Sample" uniformly, stratified sample on discrete distribution
X = zeros(nx,ny,N);
tmpVec = [-ones(1,floor(N/2)) ones(1,ceil(N/2))];
for(iX = 1:nx)
    for(iY = 1:ny)
        X(iX,iY,:)=tmpVec(randperm(N)); 
    end
end
lZ0 = nx*ny*log(2); % i0(h) should be added if h != 0

% Run sampler, adaptive or deterministic schedule
if(adaptive)
    [X,lZ,logW,alpha_log,ess_log] = smc_adpt(X,0,0,Jv,Jh,0,par);
    alpha_log = cumsum(alpha_log);
else
    [X,lZ,logW,ess_log] = smc(X,0,0,Jv,Jh,alpha,par);
    alpha_log = alpha;
end

% Compute final weights and lZ
maxlW = max( logW );
w = exp( logW - maxlW );
lZ = lZ0 + lZ + maxlW + log( sum( w ) ) - log( N );
W = w/sum(w);
end