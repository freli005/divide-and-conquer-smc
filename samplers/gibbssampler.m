function [X, Ehat] = gibbssampler(Jv, Jh, N, thinning, par)

% Blocking; update variables B x B (no overlap)
B = par.blocksize;

% Size of grid
M = size(Jv,2);

% For free boundary conditions, assume periodic with zero interaction
if(size(Jv,1) < M)
    Jv = [Jv ; zeros(1,M)];
    Jh = [Jh zeros(M,1)];
end

nx = M;
ny = M;

% Initialise
X = zeros(nx,ny,N);
Ehat = zeros(1,N);
Xnow = 2*(rand(nx,ny)<0.5)-1;

if(par.plotOn)
    figure(1);
    imagesc(Xnow,[-1 1]); colorbar;
    drawnow;
end

for(j = 1:N)
    if(par.printOn)
        fprintmod(j,100);
    end
    
    Ehat_tmp = 0;
    for(k = 1:thinning)
        % MH single flip if we have B=1
        if(B==1)
            Xnow = gibbsgrid(Xnow,Jv,Jh);
        else
            Xnow = blockgibbsgridrnd(Xnow,Jv,Jh,B,par.plotOn>1);
        end
        % Aggregate
        engy = -sum_diff_x(Xnow,Jv,Jh);
        Ehat_tmp = ((k-1)*Ehat_tmp + engy)/k;
    end
    X(:,:,j) = Xnow;
    Ehat(j) = Ehat_tmp;
    
    
    if(par.plotOn)
        figure(1);
        imagesc(Xnow,[-1 1]); colorbar;
        drawnow;
    end
end

