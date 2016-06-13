function [X,W,lZ,totMCMC,totRes,alphaStarLog] = dcsmc(Jv, Jh, N, alpha, par)
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
% N, integer: number of samples
% alpha, 1 x 2*T cell array: Annealing sequences, where n is the number of annealing
%                  steps. alpha(1) = 0, alpha(n+1) = 1

par.sampler = 'dcsmc';
printOut = par.printOn;
adaptive = numel(alpha)==0;
if(~adaptive)
    error('dcsmc: D&C-SMC has only been implemented for adaptive annealing.')
end
M = size(Jv,2);
nx = M; ny = M;
T = log2(M); % Assuming that M is a power of 2

% For analysis 
totMCMC = zeros(nx,ny,2*T+1); % Total number of MCMC iterations at each "level" of the tree
totRes = zeros(nx,ny,2*T+1); % Total number of resampling events at each "level" of the tree
alphaStarLog = zeros(nx,ny,2*T+1); % Starting inverse temperature at each "level" of the tree

% For free boundary conditions, assume periodic with zero interaction
if(size(Jv,1) < M)
    Jv = [Jv ; zeros(1,M)];
    Jh = [Jh zeros(M,1)];
end

% Keep track of the current (product) target distribution
Jv_now = zeros(nx,ny);
Jh_now = zeros(nx,ny);

% "Sample" uniformly, stratified sample on discrete distribution
X = zeros(nx,ny,N);
tmpVec = [-ones(1,floor(N/2)) ones(1,ceil(N/2))];
for(iX = 1:nx)
    for(iY = 1:ny)
        X(iX,iY,:)=tmpVec(randperm(N)); 
    end
end
logW = zeros(M,M,N);
lZ = nx*ny*log(2); % i0(h) should be added if h != 0

% Loop over the computational tree (height T)
for(k = 1:T)
    if(printOut)
        fprintf('\n\n%i\n',k);
    end
    periodic = (k==T);
    
    % Number of active population (in each direction) _after_ merge
    nActivePop = 2^(T-k);
    logWnew = zeros(nActivePop,nActivePop,N);
    % Size of populations...
    sz0 = 2^(k-1); % ...before merge
    sz1 = 2^k; % ...after merge
    
    % Masks, left/right and top/bottom, vertical and horizontal
    Mv_lr = edge_mask( sz0, sz1, par.mask_order, true, periodic, false);
    Mh_lr = edge_mask( sz0, sz1, par.mask_order, false, periodic, false);
    Mv_tb = edge_mask( sz1, sz1, par.mask_order, true, periodic, periodic);
    Mh_tb = edge_mask( sz1, sz1, par.mask_order, false, periodic, periodic);
    
    % Loop over all the subpopulations
    for(iX = 1:nActivePop)
        % row indices of subpopulations
        IX1 = ( (iX-1)*2^k+1 ) : ( (2*iX-1)*2^(k-1) );
        IX2 = ( (2*iX-1)*2^(k-1)+1 ) : ( iX*2^k );
        IX = [IX1 IX2];

        for(iY = 1:nActivePop)
            % column indices of subpopulations
            IY1 = ( (iY-1)*2^k+1 ) : ( (2*iY-1)*2^(k-1) );
            IY2 = ( (2*iY-1)*2^(k-1)+1 ) : ( iY*2^k );
            IY = [IY1 IY2];

            % Merge subpopulations left/right -----------------------------            
            
            % Get the initial/final interaction strengths
            if(printOut)
                fprintf('\n(%i,%i)->(%i,%i): ',2*iX-1,2*iY-1,2*iX-1,2*iY);
            end
            
            Jv0_top = Jv_now(IX1, IY);
            Jv1_top = Jv(IX1, IY).*Mv_lr;
            Jh0_top = Jh_now(IX1, IY);
            Jh1_top = Jh(IX1, IY).*Mh_lr;
            % Merge step / anneal (the input to the annealing is always an unweighted particle system)
            [Xtop, lZtop_merge, alpha_star_top] = smc_merge(X(IX1,IY1,:), X(IX1,IY2,:), logW(2*iX-1,2*iY-1,:), logW(2*iX-1,2*iY,:), Jv0_top, Jh0_top, Jv1_top, Jh1_top, true, par);
            [Xtop, lZtop_anneal, logWtop, alpha_log,ess_log] = smc_adpt(Xtop, Jv0_top, Jh0_top, Jv1_top, Jh1_top, alpha_star_top, par);
            lZ = lZ + lZtop_merge + lZtop_anneal;
            % Update tracking of the target distribution
            Jv_now(IX1,IY) = Jv1_top;
            Jh_now(IX1,IY) = Jh1_top;
                                    
            if(printOut)
                fprintf('%1.4f ',alpha_log);
            end
            
            % For post analysis
            totMCMC(IX1,IY,2*k) = (length(alpha_log)-2);
            totRes(IX1,IY,2*k) = sum(ess_log(1:end-1) <= par.ESS_threshold);
            alphaStarLog(IX1,IY,2*k) = alpha_star_top;
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Get the initial/final interaction strengths           
            if(printOut)
                fprintf('\n(%i,%i)->(%i,%i): ',2*iX,2*iY-1,2*iX,2*iY);
            end
            
            Jv0_bottom = Jv_now(IX2, IY);
            Jv1_bottom = Jv(IX2, IY).*Mv_lr;
            Jh0_bottom = Jh_now(IX2, IY);
            Jh1_bottom = Jh(IX2, IY).*Mh_lr;            
            % Merge step / anneal (the input to the annealing is always an unweighted particle system)
            [Xbottom, lZbottom_merge, alpha_star_bottom] = smc_merge(X(IX2,IY1,:), X(IX2,IY2,:), logW(2*iX,2*iY-1,:), logW(2*iX,2*iY,:), Jv0_bottom, Jh0_bottom, Jv1_bottom, Jh1_bottom, true, par);
            [Xbottom, lZbottom_anneal, logWbottom,alpha_log,ess_log] = smc_adpt(Xbottom, Jv0_bottom, Jh0_bottom, Jv1_bottom, Jh1_bottom, alpha_star_bottom, par);
            lZ = lZ + lZbottom_merge + lZbottom_anneal;
            % Update tracking of the target distribution
            Jv_now(IX2,IY) = Jv1_bottom;
            Jh_now(IX2,IY) = Jh1_bottom;
            
            if(printOut)
                fprintf('%1.4f ',alpha_log);
            end
                        
            % For post analysis
            totMCMC(IX2,IY,2*k) = (length(alpha_log)-2);
            totRes(IX2,IY,2*k) = sum(ess_log(1:end-1) <= par.ESS_threshold);
            alphaStarLog(IX2,IY,2*k) = alpha_star_bottom;

            
            % Merge top/bottom --------------------------------------------
            if(printOut)
                fprintf('\n(%i,%i:%i)->(%i,%i:%i): ',2*iX-1,2*iY-1,2*iY,2*iX,2*iY-1,2*iY);
            end            
            % Get the initial/final interaction strengths
            Jv0_tmp = Jv_now(IX, IY);
            Jv1_tmp = Jv(IX, IY).*Mv_tb;
            Jh0_tmp = Jh_now(IX, IY);
            Jh1_tmp = Jh(IX, IY).*Mh_tb;
            % Merge step / anneal (the input to the annealing is always an unweighted particle system)
            [Xtmp, lZtmp_merge, alpha_star_tmp] = smc_merge(Xtop, Xbottom, logWtop, logWbottom, Jv0_tmp, Jh0_tmp, Jv1_tmp, Jh1_tmp, false, par);
            [Xtmp, lZtmp_anneal, logWtmp,alpha_log,ess_log] = smc_adpt(Xtmp, Jv0_tmp, Jh0_tmp, Jv1_tmp, Jh1_tmp, alpha_star_tmp, par);
            lZ = lZ + lZtmp_merge + lZtmp_anneal;
            % Update tracking of the target distribution
            Jv_now(IX,IY) = Jv1_tmp;
            Jh_now(IX,IY) = Jh1_tmp;

            if(printOut)
                fprintf('%1.4f ',alpha_log);
            end

            % For post analysis
            totMCMC(IX,IY,2*k+1) = (length(alpha_log)-2);
            totRes(IX,IY,2*k+1) = sum(ess_log(1:end-1) <= par.ESS_threshold);
            alphaStarLog(IX,IY,2*k+1) = alpha_star_tmp;
                        
            % Update the grids
            X(IX, IY, :) = Xtmp;
            logWnew(iX, iY, :) = logWtmp;                                    
        end
    end
    
    logW = logWnew;
end

% Final weighting
maxlW = max( logW );
w = exp( logW - maxlW );
lZ = lZ + maxlW + log( sum( w ) ) - log( N );
W = w/sum(w);

% The way in which totMCMC is computed can result in negativ values, but
% these should instead be 0
totMCMC = max(totMCMC,0);
end
