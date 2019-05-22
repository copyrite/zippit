function w = uhp2udisk( z, varargin )
%UHP2UDISK Upper half-plane to unit disk
%   Map each entry of the complex matrix z by taking the upper half-plane
%   conformally to the unit disk.
%
%   uhp2udisk(z, toZero) maps z as defined above, and takes the optional
%   argument 'toZero' to origin. If omitted, toZero defaults to 1i.

p = inputParser;
addRequired(p, 'z', @(x) (isnumeric(x)));
addOptional(p, 'toZero', 1i, @(x) (numel(x) == 1 && isnumeric(x) && (imag(x) > 0)));
parse(p, z, varargin{:});

w = (p.Results.z - p.Results.toZero) ./ (p.Results.z - conj(p.Results.toZero));

end