function z = udisk2uhp( w, varargin )
%UDISK2UHP Unit disk to upper half-plane
%   Map each entry of the complex matrix z by taking the unit disk
%   conformally to the upper half-plane.
%
%   udisk2uhp(z, fromZero) maps z as defined above, and takes origin to
%   the optional argument 'fromZero'. If omitted, fromZero defaults to 1i.

p = inputParser;
addRequired(p, 'w', @(x) (isnumeric(x)));
addOptional(p, 'fromZero', 1i, @(x) (numel(x) == 1 && isnumeric(x) && (imag(x) > 0)));
parse(p, w, varargin{:});

z = (p.Results.w * conj(p.Results.fromZero) - p.Results.fromZero) ./ (p.Results.w - 1);

end