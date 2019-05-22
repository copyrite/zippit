function mz = zippit( varargin )
%ZIPPIT Streamlined constructor of a Zippit object
%   ZIPPIT() returns an empty Zippit instance.
%
%   ZIPPIT(polygon) solves the parameters required to perform a Zipper-like
%   algorithm on the given polygon, and returns a Zippit instance equipped
%   with functions that can be invoked to compute forward and inverse.
%
%   'polygon' is a complex vector of the vertices of a Jordan polygon in
%   positive orientation (i.e. counterclockwise corners of a
%   non-self-intersecting polygon).
%
%   When invoked like this, the function fills in the fields of a Zippit
%   instance. Most importantly:
%   'forward' is an approximation of a conformal function that maps the
%   polygon's interior onto the upper half-plane.
%   'inverse' is an approximation of the inverse function of 'forward'.
%   (See @Zippit/Zippit.m for details.)
%
%   ZIPPIT(polygon, algorithm) solves the parameters required to perform
%   the specified Zipper-like algorithm on the given polygon, and returns
%   a Zippit instance equipped with functions that can be invoked to
%   compute forward and inverse.
%
%   'polygon', 'forward' and 'inverse' are as defined above.
%   'algorithm' is one of the Zipper-like algorithms:
%
%       'Geodesic'
%       'Slit'
%       'Zipper'
%
%   These strings may be truncated; e.g. 'Geod' and 'Geodesic' are interpreted
%   as the same algorithm specifications.
%
%   For a detailed description of each algorithm, see:
%   "Convergence of the Zipper algorithm for conformal mapping"
%   Donald E. Marshall & Steffen Rohde
%   http://arxiv.org/abs/math/0605532.

    mz = Zippit();

    p = inputParser;
    addOptional(p, 'polygon', [], @(x) (isnumeric(x) && (numel(x) > 2)));
    addOptional(p, 'algorithm', 'Geodesic', @(x) any(validatestring(x, {'Geodesic', 'Slit', 'Zipper'})));
    parse(p, varargin{:});

    if (isempty(p.Results.polygon))
        return;
    end

    mz.polygon = reshape(p.Results.polygon, [numel(p.Results.polygon) 1]);
    mz.nv = numel(mz.polygon);
    mz.zn = nan(mz.nv+1, 1);

    if (startsWith('Geodesic', p.Results.algorithm))
        mz.geodesic();
    elseif (startsWith('Slit', p.Results.algorithm))
        mz.slit();
    elseif (startsWith('Zipper', p.Results.algorithm))
        mz.zipper();
    end

end
