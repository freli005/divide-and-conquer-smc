function X = gibbsgrid(X, Jv, Jh)
% One iteration of a single-site, random scan Gibbs sampler for the
% Ising model with periodic boundary condition.
%
% Xnow, M x M: Current state of the Markov chain
% Jv, M x M: vertical interactions, Jv(i,j) is the interaction from X(i,j)
%            to X(i+1,j). The M'th row of Jv containts the boundary
%            interactions between X(M,j) and X(1,j).
% Jh, M x M: horizontal interactions, Jh(i,j) is the interaction from
%            X(i,j) to X(i,j+1). The M'th column of Jh contains the
%            boundary interactions.
%
% For use within AIS/SMC samplers, the function supports input argument X
% to be of dimension M x M x N, in which case N updates are made in
% parallel.


[nx, ny, K] = size(X);

% ngbX: K x 4 are the 4 neighbours
% ngbJ: 1 x 4 are the corresponding interactions
indices = randperm(nx*ny);
%indices = 1:(nx*ny);

for(ii = 1:nx*ny)
    [iX, iY] = ind2sub([nx, ny], indices(ii));

    ngbX = zeros(K,4);
    ngbJ = zeros(1,4);
    % Pick out the neighbouring X:s and the corresponding edges
    if(iX == 1)
        ngbX(:,1) = X(nx,iY,:);
        ngbJ(1) = Jv(nx,iY);
    else
        ngbX(:,1) = X(iX-1,iY,:);
        ngbJ(1) = Jv(iX-1,iY);
    end
    
    if(iX == nx)
        ngbX(:,2) = X(1,iY,:);
    else
        ngbX(:,2) = X(iX+1,iY,:);
    end
    ngbJ(2) = Jv(iX,iY);
    
    if(iY == 1)
        ngbX(:,3) = X(iX,ny,:);
        ngbJ(3) = Jh(iX,ny);
    else
        ngbX(:,3) = X(iX,iY-1,:);
        ngbJ(3) = Jh(iX,iY-1);
    end
    
    if(iY == ny)
        ngbX(:,4) = X(iX,1,:);
    else
        ngbX(:,4) = X(iX,iY+1,:);
    end
    ngbJ(4) = Jh(iX, iY);
    
    % Fix double counting for 1d models
    if(nx == 1), ngbJ(2) = 0; end;
    if(ny == 1), ngbJ(4) = 0; end;
    
    % Gibbs sampling
%     logPp = sum(bsxfun(@times, ngbJ, ngbX),2);
%     logPm = sum(bsxfun(@times, ngbJ, -ngbX),2);    
%     C = max(logPp, logPm);
%     Pm = exp(logPm-C)./(exp(logPm-C)+exp(logPp-C)); % P(X = -1)
%     X(iX,iY,:) = 1 - 2*(rand(1,1,K) < Pm);
    
    % Flip HM
    Xnow = X(iX,iY,:);
    logAP = sum(bsxfun(@times, ngbJ, bsxfun(@times, ngbX, -Xnow(:))),2) - sum(bsxfun(@times, ngbJ, bsxfun(@times, ngbX, Xnow(:))),2);
    ind = rand(K,1) < exp(logAP);
    if(any(ind))
        X(iX,iY,ind) = -Xnow(ind);
    end  
end