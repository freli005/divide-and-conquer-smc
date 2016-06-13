function X = blockgibbsgridrnd(X, Jv, Jh, B, plotOn)
% One iteration of a Blocked, random scan Gibbs sampler for the
% Ising model with periodic boundary condition.
%
% Xnow, M x M: Current state of the Markov chain
% Jv, M x M: vertical interactions, Jv(i,j) is the interaction from X(i,j)
%            to X(i+1,j). The M'th row of Jv containts the boundary
%            interactions between X(M,j) and X(1,j).
% Jh, M x M: horizontal interactions, Jh(i,j) is the interaction from
%            X(i,j) to X(i,j+1). The M'th column of Jh contains the
%            boundary interactions.
% B: Block size
%
% For use within AIS/SMC samplers, the function supports input argument X
% to be of dimension M x M x N, in which case N updates are made in
% parallel.

[M, ny, K] = size(X);

if(M~=ny)
    error('Grid has to be square');
end

m = M/B; % Number of blocks in each direction (on average)

if(m <= 1)
    error('Need at least 2 blocks');
end

numIter = round(m^2); % Run this many iterations of a random scan (i.e. not random permutation) Gibbs sampler

% Trial blocks used in Gibbs sampler
numPossibilities = 2^(B^2);
if(numPossibilities > 1e6)
    error('Are you sure...?');
end
trialBlocks = combn([-1 1], B^2)'; % [B^2, numPosibilities]
trialBlocks = reshape(trialBlocks, [B B 1 numPossibilities]); % Put numPossibilities on fourth dimension since K is on third
% Precompute, to use for within block interactions
horzProd = trialBlocks(:,1:end-1,:,:).*trialBlocks(:,2:end,:,:);
vertProd = trialBlocks(1:end-1,:,:,:).*trialBlocks(2:end,:,:,:);
% Precomupte, to use for between block interactions
edgeX = zeros(B,4,1,numPossibilities);
% top bottom left right
edgeX(:,1,1,:) = trialBlocks(1,:,1,:);
edgeX(:,2,1,:) = trialBlocks(end,:,1,:);
edgeX(:,3,1,:) = trialBlocks(:,1,1,:);
edgeX(:,4,1,:) = trialBlocks(:,end,1,:);

% ngbX: B x 4 x K are the 4*B neighbours
% ngbJ: B x 4 are the corresponding interactions
indices = randi(M^2, [numIter 1]);

for(ii = 1:numIter)
    % Indices to the top left corner of the block (chosen randomly)
    [iX, iY] = ind2sub([M, M], indices(ii));

    % Compute indices (rows and cols) to the current block
    iXX = mod( (iX:(iX+B-1))-1, M)+1;
    iYY = mod( (iY:(iY+B-1))-1, M)+1;
   
    ngbX = zeros(B,4,K);
    ngbJ = zeros(B,4);
    
    % Pick out the neighbouring X:s and the corresponding edges
    if(iX == 1) % top border (wrapped)
        ngbX(:,1,:) = X(end, iYY, :);
        ngbJ(:,1) = Jv(end, iYY);
    else % top border
        ngbX(:,1,:) = X(iX-1, iYY, :);
        ngbJ(:,1) = Jv(iX-1, iYY);
    end
    
    if(iXX(end) == M) % bottom (wrapped)
        ngbX(:,2,:) = X(1, iYY, :);
    else % bottom
        ngbX(:,2,:) = X(iXX(end)+1, iYY, :);
    end
    ngbJ(:,2) = Jv(iXX(end), iYY);
    
    if(iY == 1) % left (wrapped)
        ngbX(:,3,:) = X(iXX, end, :);
        ngbJ(:, 3) = Jh(iXX, end);
    else % left
        ngbX(:,3,:) = X(iXX, iY-1, :);
        ngbJ(:,3) = Jh(iXX, iY-1);
    end
    
    if(iYY(end) == M) % right (wrapped)
        ngbX(:,4,:) = X(iXX, 1, :);
    else % right
        ngbX(:,4,:) = X(iXX, iYY(end)+1, :);
    end
    ngbJ(:,4) = Jh(iXX, iYY(end));
    
    %%% Gibbs sampling
    
    % Compute within-block-interactions
    logP = sum(sum(bsxfun(@times, Jh(iXX, iYY(1:end-1)), horzProd),1),2) + ...
        sum(sum(bsxfun(@times, Jv(iXX(1:end-1), iYY), vertProd),1),2);
    
    % Compute neighbour-interactions
    logP = logP + sum(sum(bsxfun(@times, ngbJ, bsxfun(@times, ngbX, edgeX)),1),2);
    C = max(logP,[], 4);
    logP = bsxfun(@minus, logP, C);
    bins = cumsum( shiftdim(exp(logP),2), 2 ); % [1 1 K numPossibilities]  -->  [K numPossibilities]
    bins = bsxfun(@rdivide, bins, bins(:,end)); 
    inds = numPossibilities + 1 - sum(bsxfun(@lt, rand(K,1), bins),2);
    
    X(iXX,iYY,:) = trialBlocks(:,:,inds);
 
    if(plotOn && K==1)
        imagesc(X);
        title(ii);
        drawnow;
    end
end