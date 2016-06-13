function [logW, alpha, ess, ITER] = adpt_alpha(S, logW, alpha_init, alpha_max, cess_target, epsilon)
% SMC utility function - adapts the value of the annealing parameter ALPHA
% based on the CESS criterion.
% Input:
%   S : [N,1] Negative energy of (new) target distribution, 
%       \pi(x^i) \propto exp(alpha*S^i), for each particle x^i.
%   logW : [N,1] log-importance weights of current sample.
%   alpha_init : initialise search procedure at this value.
%   alpha_max : search for alpha in [0, alpha_max].
%   cess_target : target CESS in interval [0,1]
%
% The function computes and returns the largest value of alpha \in
% [0,alpha_max] such that CESS >= cess_target, if such a value exists, or
% returns alpha = 0 otherwise. The function also returns the updated
% log-weights and the ESS for the found value of alpha. The function uses a
% combination newton and/or bisection method.

N = length(S);

% Subtract maximums / =divide by common constants
maxlW = max(logW);
maxS = max(S);
lWd = logW -maxlW;
Sd = S-maxS;

% Normalisation for the weights
WN = sum(exp(lWd));

% Compute CESS at alpha_max
F1 = sum(exp(lWd + alpha_max*Sd));
F2 = sum(exp(lWd + 2*alpha_max*Sd));
cess = F1^2/F2/WN;

if(cess >= cess_target)
    alpha = alpha_max;
    ITER = 0;
else
    ITER = 1;
        
    % Compute CESS at alpha_init
    alpha = alpha_init;
    F1 = sum(exp(lWd + alpha*Sd));
    F2 = sum(exp(lWd + 2*alpha*Sd));
    cess = F1^2/F2/WN;

    %%% Try first with Newton's method - at instability, switch to bisection
    newton = 1;
    while(abs(cess-cess_target) >= epsilon) % While tolerance threshold not reached
        ITER = ITER+1;
        if(ITER > 10000)
            error('SMC-Utility, adpt-alpha: Bisection method failed in 10000 iterations');
        end
        
        if(newton)
            % Compute Derivative
            Fs1 = sum(Sd.*exp(lWd + alpha*Sd));
            Fs2 = sum(Sd.*exp(lWd + 2*alpha*Sd));
            derivative = 2*F1/F2^2/WN*(Fs1*F2 - Fs2*F1);
            
            % Update alpha
            alpha = alpha - (cess-cess_target)/derivative;
            if(alpha < 0 || alpha > alpha_max)
                newton = 0;
                a = 0;
                b = alpha_max;
                alpha = (a+b)/2; % Reset at initial value (use mid-point for bisection)
            end
        else
            % Choose left or right side
            if(cess > cess_target) % Increase alpha -> continue with right interval
                a = alpha;
            else % Decrease alpha -> continue with left interval
                b = alpha;
            end
            % Bisect interval
            alpha = (a + b)/2;
        end
        
        % Compute new CESS
        F1 = sum(exp(lWd + alpha*Sd));
        F2 = sum(exp(lWd + 2*alpha*Sd));
        cess = F1^2/F2/WN;
    end
end

% Compute ESS and updated log-weights
ess = F1^2/sum(exp(2*lWd+2*alpha*Sd))/N;
logW = logW + alpha*S;
