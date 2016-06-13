function [X,lZ,alpha] = smc_merge(X1, X2, logW1, logW2, Jv0, Jh0, Jv1, Jh1, leftright, par)
% D&C-SMC utility function:
% Given two weighted set of particles, such that the product measure
% targets J(0), generate one weighted set of particles targeting J(alpha*),
%
%   J(alpha) = J(0) + alpha*J(1),
%
% where alpha* is chosen adaptively.
%
% par.merge_method    controls which method that is used. 
%   0: D&C-SIR, i.e. sample independently for each subpopulation and set
%      alpha* = 1.
%   1: D&C-SMC (mix), where alpha* is adapted based on the ESS of the N^2
%      particle system. The threshold is ESS(N^2 system) >= t*N^a where
%      t = par.ESS_merge_threshold, and a = par.merge_exponent.
%   2: D&C-SMC (mix), where alpha* is adapted based on the CESS of the 2
%      marginal systems (each comprising N particles). The threshold is
%      CESS(each marginal N system) >= t*N.
%   3: D&C-SMC (lw), where alpha* is adapted based similarly to case 1.
%
% For case (mix), the N^2 particle system is subsampled using a stratified
% resampling procedure, either directly for the N^2 particle system
% (par.merge_subsampling = 1) or based on the marginal/conditional thereof
% (par.merge_subsampling = 2). For cases (sir) and (lw) the value of
% parameter par.merge_subsampling is irrelevant. For case (sir) we always
% resample each subpopulation independently. For case (lw) we subsample the
% mN particle system with stratified resampling (corresponding to
% par.merge_subsampling = 1).
%
% The function uses a newton and/or bisection method to find alpha*.

if(par.merge_method == 0) % SIR style
    [X,lZ,alpha] = smc_merge_sir(X1, X2, logW1, logW2, leftright);
elseif(par.merge_method == 1 || par.merge_method == 2) % (mix)
    [X,lZ,alpha] = smc_merge_mix(X1, X2, logW1, logW2, Jv0, Jh0, Jv1, Jh1, leftright, par);
elseif(par.merge_method == 3); % (lw)
    [X,lZ,alpha] = smc_merge_lw(X1, X2, logW1, logW2, Jv0, Jh0, Jv1, Jh1, leftright, par);
else
    error('smc_merge: unknown merge method');
end