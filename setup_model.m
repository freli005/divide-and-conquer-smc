function [Jv,Jh] = setup_model(n)
% This function defines the parameters of the Ising model. It is called from
% the various sampling routines.

periodic = true;                    % Periodic boundary condition (yes/no)
M = 2^n;                            % Size of square lattice (MxM) - needs to be a multiple of 2
inverseTemperature = 1/2.269185;    % Set the inverse temperature for the model

% Construct matrices of edge potentials. Jv contains all vertical edges and
% Jh contains all horizontal edges. Note that it is possible to use an
% inhomogeneous temperature if desired, by appropriate changes to the
% defintions of Jv and Jh.
if(periodic)
    Jv = inverseTemperature*ones(M);
    Jh = inverseTemperature*ones(M);
else
    Jv = inverseTemperature*ones(M-1,M);
    Jh = inverseTemperature*ones(M,M-1);
end
