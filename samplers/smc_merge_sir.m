function [X,lZ,alpha] = smc_merge_sir(X1, X2, logW1, logW2, leftright)
% D&C-SMC utility function:
% Given two weighted set of particles, such that the product measure
% targets J(0), generate one weighted set of particles targeting J(0),
% by independently resampling each subpopulation.

[nx,ny,N] = size(X1);

% normalise the weights for the two marginal particle systems
logW1 = logW1(:);
maxlW1 = max( logW1 );
lWd1 = logW1 - maxlW1;
w1 = exp(lWd1);
C1 = sum(w1);
W1 = w1/C1;

%
logW2 = logW2(:)';
maxlW2 = max( logW2 );
lWd2 = logW2 - maxlW2;
w2 = exp(lWd2);
C2 = sum(w2);
W2 = w2/C2;

% Resample each system independently
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

