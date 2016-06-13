function fprintmod(k,n)
% FPRINTMOD(K,N) prints K to prompt if K divides N.

if(~mod(k,n))
    fprintf('%i\n',k);
end