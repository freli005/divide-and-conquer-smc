function [I,J] = stratified_resampling_N2toN(Wij,N)
% Wij : [N N], normalised weights

Wti = sum(Wij,2); % Marginal row probabilities

% Stratified sampling of rows
qc = cumsum(Wti);
u = ((0:N-1)+rand(1,N))/N;

[~,ind1]=sort([u(:) ; qc(:)]);
ind2=find(ind1<=N);
I = ind2'-(0:N-1);

% Conditional stratified sampling of columns
qc = bsxfun(@rdivide, cumsum(Wij,2), Wti); % CDF for every row
u = ((0:N-1)+rand(1,N))/N;

sigma = randperm(N);
J = zeros(N,1);
for(i = 1:N)
    [~,ind1]=sort([u(sigma(i)) qc(I(i),:)]);
    ind2=find(ind1==1);
    J(i) = ind2;
end

end