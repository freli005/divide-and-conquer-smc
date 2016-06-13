function [X,lZ,alpha] = smc_merge_mix(X1, X2, logW1, logW2, Jv0, Jh0, Jv1, Jh1, leftright, par)
% D&C-SMC utility function:
% Given two weighted set of particles, such that the product measure
% targets J(0), generate one weighted set of particles targeting J(alpha*),
%
%   J(alpha) = J(0) + alpha*J(1),
%
% where alpha* is chosen adaptively.
%
% par.merge_method    controls which method that is used. 
%   1: D&C-SMC (mix), where alpha* is adapted based on the ESS of the N^2
%      particle system. The threshold is ESS(N^2 system) >= t*N^a where
%      t = par.ESS_merge_threshold, and a = par.merge_exponent.
%   2: D&C-SMC (mix), where alpha* is adapted based on the CESS of the 2
%
% For case (mix), the N^2 particle system is subsampled using a stratified
% resampling procedure, either directly for the N^2 particle system
% (par.merge_subsampling = 1) or based on the marginal/conditional thereof
% (par.merge_subsampling = 2).
%
% The function uses a newton and/or bisection method to find alpha*.

[nx,ny,N] = size(X1);
ess_target = par.ESS_merge_threshold; % in [0,1], or larger if exponent < 2

% In the CESS case, make sure that the threshold is in the correct interval
if(ess_target < 0 || (par.merge_method == 2 && ess_target > 1))
    error('(C)ESS target incorrectly set');
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
ess1 = 1/sum(W1.^2); % ESS in [1,N]
%
logW2 = logW2(:)';
maxlW2 = max( logW2 );
lWd2 = logW2 - maxlW2;
w2 = exp(lWd2);
C2 = sum(w2);
W2 = w2/C2;
ess2 = 1/sum(W2.^2); % ESS in [1,N]

% If it is not possible to satisfy the (C)ESS target we get alpha=0, which
% means that the mixture sampling method reduces to the SIR style sampling.
if(par.merge_method == 1 && ess1*ess2 <= ess_target*N^(par.merge_exponent))
    ind1 = resampling(W1,3); % Stratified!
    ind2 = resampling(W2,3); % Stratified!
    ind1 = ind1(randperm(N));
    ind2 = ind2(randperm(N));
    
    % Output
    lZ = log(sum(w1)) + log(sum(w2)) + maxlW1 + maxlW2 - 2*log(N);
    alpha = 0;
        
    if(leftright)
        X = [X1(:,:,ind1) X2(:,:,ind2)];
    else
        X = [X1(:,:,ind1) ; X2(:,:,ind2)];
    end
    
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expand N^2 particles and compute S
if(par.mask_order ~= 1)
    % For masko_order == 1 a much faster implementation is use below
    S = zeros(N);
    if(leftright) % left/right
        for(i = 1:N)
            XX = [repmat( X1(:,:,i),[1 1 N]) X2];
            S(i,:) = sum_diff_x(XX,Jvt,Jht);
        end
    else
        for(i = 1:N)
            XX = [repmat( X1(:,:,i),[1 1 N]) ; X2];
            S(i,:) = sum_diff_x(XX,Jvt,Jht);
        end
    end
    %SS = S;
else
    % This is a hard coded version of the above for a specific structure;
    % only works when: "gibbsmask" is all ones; Jh,Jv are homogeneous (can
    % easily be generalise); we have the specific energy function of the Ising
    % model
    
    S = zeros(N);
    if(leftright)
        % Only Jht will have non-zero entries
        %       [0 0 ... const 0 ... 0]           [0 ... 0 const 0 ... 0 const]
        % Jht = [0 0 ... const 0 ... 0]  or Jht = [0 ... 0 const 0 ... 0 const]
        %       [0 0 ... const 0 ... 0]           [0 ... 0 const 0 ... 0 const]
        bp = find(Jht(1,:));
        periodic = numel(bp)>1;
        for(k = 1:nx)
            Xd = squeeze( X1(k,end,:) ) *squeeze( X2(k,1,:) )';
            S = S + Jht(k,bp(1))*Xd;
            if(periodic)
                Xd = squeeze( X1(k,1,:) )*squeeze( X2(k,end,:) )';
                S = S + Jht(k,bp(2))*Xd;
            end
        end
    else
        % Only Jvt will have non-zero entries
        bp = find(Jvt(:,1));
        periodic = numel(bp)>1;
        for(k = 1:ny)
            Xd = squeeze(X1(end,k,:)) * squeeze(X2(1,k,:))';
            S = S + Jvt(bp(1),k)*Xd;
            if(periodic)
                Xd = squeeze(X1(1,k,:)) * squeeze(X2(end,k,:))';
                S = S + Jvt(bp(2),k)*Xd;
            end
        end
    end
end

% Initial joint log-weights
logW = bsxfun(@plus, logW1, logW2); % N^2 initial log-weights

% Subtract maximums / =divide by common constants
maxlW = max(logW(:));
maxS = max(S(:));
lWd = logW-maxlW;
Sd = S-maxS;

%%%% Use newton/bisection to find alpha* %%%%
epsilon = par.CESS_tolerance;
a = 0; b = 1; % Start and end points of the interval

% Compute at alpha_max = 1
alpha = 1;
if(par.merge_method == 1) % ESS-based
    ess_target = ess_target*N^(par.merge_exponent-2); % In [0,1]
    
    num = (sum(sum(exp(lWd + alpha*Sd))))^2;
    den = sum(sum(exp(2*lWd+2*alpha*Sd)));
    essN2 = num/den/N^2; % In [0,1]
    
    if(essN2 < ess_target)
        % Start at the mid-point
        alpha = (a + b)/2;
        num = (sum(sum(exp(lWd + alpha*Sd))))^2;
        den = sum(sum(exp(2*lWd+2*alpha*Sd)));
        essN2 = num/den/N^2; % In [0,1]
        
        while(abs(essN2-ess_target) >= epsilon) % While tolerance threshold not reached
            % Choose left or right side
            if(essN2 > ess_target) % Increase alpha -> continue with right interval
                a = alpha;
            else % Decrease alpha -> continue with left interval
                b = alpha;
            end
            % Bisect interval
            alpha = (a + b)/2;
            num = (sum(sum(exp(lWd + alpha*Sd))))^2;
            den = sum(sum(exp(2*lWd+2*alpha*Sd)));
            essN2 = num/den/N^2; % In [0,1]
        end
    end
else % Marginal CESS-based
    % pre-compute some alpha-independent quantities that are used to comptue CESS
    lWd_half1 = bsxfun(@minus, lWd, lWd1/2);
    lWd_half2 = bsxfun(@minus, lWd, lWd2/2);
    
    % CESS at alpha_max
    num = (sum(sum(exp(lWd + alpha*Sd))))^2;
    den1 = sum((sum(exp(lWd_half1 + alpha*Sd),2)).^2, 1);
    den2 = sum((sum(exp(lWd_half2 + alpha*Sd),1)).^2, 2);
    cess1 = num/den1/C1; % In [0,1]
    cess2 = num/den2/C2; % In [0,1]
    
    if(cess1 < ess_target-epsilon || cess2 < ess_target-epsilon)
        % Start at the mid-point
        alpha = (a + b)/2;
        num = (sum(sum(exp(lWd + alpha*Sd))))^2;
        den1 = sum((sum(exp(lWd_half1 + alpha*Sd),2)).^2, 1);      
        den2 = sum((sum(exp(lWd_half2 + alpha*Sd),1)).^2, 2);        
        cess1 = num/den1/C1; % In [0,1]
        cess2 = num/den2/C2; % In [0,1]
        
        % Both must be bigger (up to tolerance), and at least
        % one must be within absolue tolerance
        bothAreBigger  = (cess1 > ess_target-epsilon) && (cess2 > ess_target-epsilon);
        done = bothAreBigger && (cess1 < ess_target + epsilon || cess2 < ess_target + epsilon);
        
        while(~done) % While tolerance threshold not reached            
            % Choose left or right side
            if(bothAreBigger) % Increase alpha -> continue with right interval
                a = alpha;
            else % Decrease alpha -> continue with left interval
                b = alpha;
            end
            % Bisect interval
            alpha = (a + b)/2;
            num = (sum(sum(exp(lWd + alpha*Sd))))^2;
            den1 = sum((sum(exp(lWd_half1 + alpha*Sd),2)).^2, 1);
            den2 = sum((sum(exp(lWd_half2 + alpha*Sd),1)).^2, 2);
            cess1 = num/den1/C1; % In [0,1]
            cess2 = num/den2/C2; % In [0,1]
            
            bothAreBigger  = (cess1 > ess_target-epsilon) && (cess2 > ess_target-epsilon);
            done = bothAreBigger && (cess1 < ess_target + epsilon || cess2 < ess_target + epsilon);            
        end
    end
end

% Compute weights
logW = logW + alpha*S;
maxlW = max(logW(:));
w = exp( logW - maxlW );
sumw = sum(sum(w));
W = w/sumw;
lZ = log(sumw) + maxlW - 2*log(N);

if(par.merge_subsampling == 2) % Stratified on marginal
    [ind1,ind2] = stratified_resampling_N2toN(W,N);
elseif(par.merge_subsampling == 1) % Stratified on joint
    ind = stratified_resampling_MtoN(W(:),N);
    [ind1,ind2] = ind2sub([N,N],ind);
else
    error('Unknown merge method');
end

% Output
if(leftright)
    X = [X1(:,:,ind1) X2(:,:,ind2)];
else
    X = [X1(:,:,ind1) ; X2(:,:,ind2)];
end
