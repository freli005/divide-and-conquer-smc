function S = sum_diff_x(X,Jv,Jh)
% SMC utility function 
% Computes sum_{ij} J_{ij}*x_i*x_j for every page in X

Xdv = X .* [X(2:end,:,:) ; X(1,:,:)];
Xdh = X .* [X(:,2:end,:)   X(:,1,:)];
S = sum(sum(bsxfun(@times, Jh, Xdh),1),2) + sum(sum(bsxfun(@times, Jv, Xdv),1),2);
S = S(:);
end