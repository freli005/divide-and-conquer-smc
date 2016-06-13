function [X,lZ,alpha] = smc_merge_lw(X1, X2, logW1, logW2, Jv0, Jh0, Jv1, Jh1, leftright, par)
% D&C-SMC utility function:
% Given two weighted set of particles, such that the product measure
% targets J(0), generate one weighted set of particles targeting J(alpha*),
%
%   J(alpha) = J(0) + alpha*J(1),
%
% where alpha* is chosen to attain a threshold ESS. The methods uses
% "lightweight sampling", i.e. it implements the following algorithm:
%  1. Sample m*N particles from each marginal system based on their
%     importance weights.
%  2. Pair the m*N particles together to form m*N particles targeting the
%     joint. Select alpha* so that the ESS of the m*N particle system is at
%     least t*N.
%
% Thus, the algorithm requires the following settings (for compatibility
% with the smc_merge_mix function),
%
%   par.merge_exponent = 1
%   par.ESS_merge_threshold <= m
%  
% Furthermore, the function is restricted to par.merge_subsampling = 1,
% meaning that the m*N particle system is subsampled using a stratified
% resampling procedure "directly".
%
% the choice of m is specified in par.lw_m.
%
% The function uses a newton/bisection method to find alpha*.

m = par.lw_m;
[nx,ny,N] = size(X1);
ess_target = par.ESS_merge_threshold; % in [0,1], or larger if exponent < 2
exponent = par.merge_exponent;

% Check settings
if( exponent ~= 1 || par.merge_subsampling ~= 1)
    error('smc_merge_lw: the lightweight mixture sampling procedure is incompatible with the chosen settings.');
end

% Make sure that the threshold is in the correct interval
if(ess_target < 1 || ess_target > m || round(m)~=m)
    error('smc_merge_lw: ESS target and/or m incorrectly set.');
end

Jvt = Jv1-Jv0; % Compute incremental edge potentials / J = J0 + alpha*Jt
Jht = Jh1-Jh0;

% Compute some statistics on the marginal particle systems
logW1 = logW1(:);
maxlW1 = max( logW1 );
lWd1 = logW1 - maxlW1;
w1 = exp(lWd1);
C1 = sum(w1);
W1 = w1/C1;
%ess1 = 1/sum(W1.^2); % ESS in [1,N]
%
logW2 = logW2(:)';
maxlW2 = max( logW2 );
lWd2 = logW2 - maxlW2;
w2 = exp(lWd2);
C2 = sum(w2);
W2 = w2/C2;
%ess2 = 1/sum(W2.^2); % ESS in [1,N]

% Sample w. replacement (stratified) on each system independently
% [This could likely be made a bit quicker by clever implementation]
ind1_MM = stratified_resampling_MtoN(W1,m*N);
ind2_MM = stratified_resampling_MtoN(W2,m*N);
ind2_MM = ind2_MM(randperm(m*N)); % Randomise one of the systems

% Contribution to lZ from initial (marginal) weights
lZ = log(sum(w1)) + log(sum(w2)) + maxlW1 + maxlW2 - 2*log(N);

% Expand m*N particles and compute S
if(par.mask_order ~= 1)
    error('smc_merge_lw: mask order > 1 has not been implemented for the lw procedure. Can easily be adapted from smc_merge_mix.');
else
    % This is a hard coded version of the above for a specific structure;
    % only works when: "gibbsmask" is all ones; Jh,Jv are homogeneous (can
    % easily be generalise); we have the specific energy function of the Ising
    % model
    
    if(leftright)
        % Only Jht will have non-zero entries
        %       [0 0 ... const 0 ... 0]           [0 ... 0 const 0 ... 0 const]
        % Jht = [0 0 ... const 0 ... 0]  or Jht = [0 ... 0 const 0 ... 0 const]
        %       [0 0 ... const 0 ... 0]           [0 ... 0 const 0 ... 0 const]
        bp = find(Jht(1,:));
        periodic = numel(bp)>1;
       
        Xd = bsxfun(@times, Jht(:,bp(1)),  X1(:,end,ind1_MM).*X2(:,1,ind2_MM) ); % [nx 1 (m*N)]
        S = squeeze(sum(Xd,1));
        if(periodic)
            Xd = bsxfun(@times, Jht(:,bp(2)),  X1(:,1,ind1_MM).*X2(:,end,ind2_MM) ); % [nx 1 (m*N)]
            S = S + squeeze(sum(Xd,1));
        end
    else
        % Only Jvt will have non-zero entries
        bp = find(Jvt(:,1));
        periodic = numel(bp)>1;
        
        Xd = bsxfun(@times, Jvt(bp(1),:),  X1(end,:,ind1_MM).*X2(1,:,ind2_MM) ); % [1 ny (m*N)]
        S = squeeze(sum(Xd,2));
        if(periodic)
            Xd = bsxfun(@times, Jvt(bp(2),:),  X1(1,:,ind1_MM).*X2(end,:,ind2_MM) ); % [1 ny (m*N)]
            S = S + squeeze(sum(Xd,2));
        end
    end
end

% Subtract maximums / =divide by common constants
maxS = max(S);
Sd = S-maxS;

%%%% Use bisection to find alpha* %%%%%%%
epsilon = par.CESS_tolerance;
a = 0; b = 1; % Start and end points of the interval

% Compute at alpha_max = 1
alpha = 1;
ess_target = ess_target/m; % In [0,1]

num = (sum(exp(alpha*Sd)))^2;
den = sum(exp(2*alpha*Sd));
essMN = num/den/(m*N); % In [0,1]

if(essMN < ess_target)
    % Start at the mid-point
    alpha = (a + b)/2;
    num = (sum(exp(alpha*Sd)))^2;
    den = sum(exp(2*alpha*Sd));
    essMN = num/den/(m*N); % In [0,1]
    
    while(abs(essMN-ess_target) >= epsilon) % While tolerance threshold not reached
        % Choose left or right side
        if(essMN > ess_target) % Increase alpha -> continue with right interval
            a = alpha;
        else % Decrease alpha -> continue with left interval
            b = alpha;
        end
        % Bisect interval
        alpha = (a + b)/2;
        num = (sum(exp(alpha*Sd)))^2;
        den = sum(exp(2*alpha*Sd));
        essMN = num/den/(m*N); % In [0,1]
    end
end

% Compute weights
logW = alpha*S;
maxlW = max(logW);
w = exp( logW - maxlW );
sumw = sum(sum(w));
W = w/sumw;
lZ = lZ + log(sumw) + maxlW - log(m*N); 

ind = stratified_resampling_MtoN(W(:),N);
ind1 = ind1_MM(ind);
ind2 = ind2_MM(ind);

% Output
if(leftright)
    X = [X1(:,:,ind1) X2(:,:,ind2)];
else
    X = [X1(:,:,ind1) ; X2(:,:,ind2)];
end
