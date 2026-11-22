function tf = isPositiveIntegerArray(x)
% ISPOSITIVEINTEGERARRAY True if x is a non-empty numeric array whose
% elements are all positive integers (e.g., 3 or [2 1 3]).
%
% Replaces the undocumented MATLAB-internal function
% isPositiveIntegerValuedNumeric, which is not guaranteed to exist in every
% MATLAB release.
    tf = isnumeric(x) && ~isempty(x) && all(isfinite(x(:))) && ...
         all(x(:) > 0) && all(x(:) == round(x(:)));
end
