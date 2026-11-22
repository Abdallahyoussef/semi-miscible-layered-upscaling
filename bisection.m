function [x, fval, exitflag] = bisection(f, lb, ub, target, options)
% BISECTION  Scalar root finding by interval bisection: f(x) = target.
%
%   [x, fval, exitflag] = bisection(f, lb, ub)
%   [x, fval, exitflag] = bisection(f, lb, ub, target)
%   [x, fval, exitflag] = bisection(f, lb, ub, target, options)
%
%   INPUTS
%     f       - function handle of one scalar variable
%     lb, ub  - lower and upper bounds that bracket the root
%     target  - value to solve for (default 0)
%     options - structure with optional fields
%                 .TolX   - stop when half the bracket width <= TolX (default 1e-12)
%                 .TolFun - stop when |f(x) - target| <= TolFun   (default 0)
%                 .MaxIter- maximum iterations                     (default 200)
%
%   OUTPUTS
%     x        - root estimate
%     fval     - f(x) - target at the returned x
%     exitflag -  1  converged on TolX
%                 2  converged on TolFun
%                 0  maximum iterations reached
%                -2  the bounds do not bracket a root
%
%   Same call interface as the MATLAB File Exchange function "Bisection
%   method root finding" (File Exchange #28150), which was used during
%   development. This file is an independent implementation so that the
%   repository has no external dependency.

    if nargin < 4 || isempty(target), target = 0; end
    if nargin < 5 || isempty(options), options = struct(); end
    TolX    = getOpt(options, 'TolX', 1e-12);
    TolFun  = getOpt(options, 'TolFun', 0);
    MaxIter = getOpt(options, 'MaxIter', 200);

    g  = @(s) f(s) - target;
    a  = lb;  b  = ub;
    fa = g(a); fb = g(b);

    if fa == 0, x = a; fval = fa; exitflag = 2; return; end
    if fb == 0, x = b; fval = fb; exitflag = 2; return; end
    if sign(fa) == sign(fb)
        x = NaN; fval = NaN; exitflag = -2;
        return;
    end

    exitflag = 0;
    for k = 1:MaxIter
        x    = (a + b) / 2;
        fval = g(x);
        if abs(fval) <= TolFun
            exitflag = 2; return;
        end
        if sign(fval) == sign(fa)
            a = x; fa = fval;
        else
            b = x;
        end
        if (b - a) / 2 <= TolX
            exitflag = 1; return;
        end
    end
end

function v = getOpt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
