function m = bisection(f,a,b,epsilon)
% BISECTION of decreasing function
%   ROOT = BISECTION(F,A,B,EPSILON) finds a root of the monotone
%   continuous function F on the interval [A,B] up to precision EPSILON.
%
%   If F does not have a root in [A,B], the value B is returned.
%
% Fredrik Lindsten
% 2015-04-25

fa = f(a);
fb = f(b);

if(fa*fb > 0) % No root, always output right value
    m = b;
else
    % Start at the mid-point
    m = (a + b)/2;
    fm = f(m);
    
    while(abs(fm) >= epsilon) % While tolerance threshold not reached
        % Choose left or right side
        if(fa*fm > 0) % Same sign to the left, continue with right interval
            a = m;
            fa = fm;
        else % Different sign to the left, continue with left interval
            b = m;
        end
        % Bisect interval
        m = (a + b)/2;
        fm = f(m);
    end
end
