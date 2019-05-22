function phi = phi_slit_fwd_i( w, d )
%PHI_SLIT_FWD_I Forward slit phi_i
%   A component function of a Zipper-like algorithm.
%
%   The forward phi_i of slit algorithm is the inverse of slit map. However,
%   the slit map has no closed form inverse. Therefore we resort to
%   applying Newton's method.
%
%   Slit map g and its derivative g' are
%
%   g (z) =   (z-p)^ p    (z+1-p)^(1-p)
%   
%   g'(z) = p (z-p)^(p-1) (z+1-p)^(1-p)  +  (z-p)^p (1-p)(z+1-p)^-p
%         =   (z-p)^(p-1) (z+1-p)^-p [ p (z+1-p) + (1-p) (z-p) ]
%         =   (z-p)^(p-1) (z+1-p)^-p [ z ] .
%   
%   Solve g (z) = w for z by Newton's method:
%   z_{n+1} = z_{n} - (g (z_{n}) - w) / g'(z_{n}) .
%   
%   Iteration is improved further by dividing the domain into four
%   subsets and in each applying a suitable "preliminary map" (see first
%   reference). Iterates that erroneously escape the upper half-plane
%   are forced back by taking the complex conjugate.
%   
%   References:
%   "Convergence of the Zipper algorithm for conformal mapping"
%   Donald E. Marshall & Steffen Rohde
%   http://arxiv.org/abs/math/0605532
%   
%   Package "ConfomalMaps" for the programming language "Julia"
%   Samuel S. Watson
%   https://github.com/sswatson/ConformalMaps.jl

    if (numel(d) ~= 1 || isnan(d) || isinf(d))
        error('Zipper error: Invalid parameter in forward slit map');
    end

    p = angle(d) / pi;
    C = abs(d) / (p^p * (1-p)^(1-p));
    
    f  = @(z) (z-p).^p .* (z+1-p).^(1-p);
    df = @(z) (z-p).^(p-1) .* (z+1-p).^-p .* z;

    % TODO adjustable tolerance
    tol = 10^-14;

    function zn = errcorr(zn)
        % Simple error-correcting function
        err = imag(zn) < 0;
        zn(err) = conj(zn(err));
    end

    function z = newton(ff, dff, z, w, corrected)
        % Keep track of points to operate upon. Immediately eliminate
        % perfectly guessed points.
        op = ((imag(z) < 0) | abs(ff(z, w)) > tol);
        
        % Whether we wish to force points onto the upper half-plane.
        % Pick and apply the desired correcting function.
        if corrected
            corrfn = @errcorr;
        else
            corrfn = @(x) x;
        end


        for ii=1:100
            % Apply Newton and correction
            z(op) = corrfn( z(op) - ff(z(op), w(op)) ./ dff(z(op), w(op)) );

            % Eliminate converged points
            op(op) = abs(ff(z(op), w(op))) > tol;

            % If everything has converged, quit early
            if isempty(find(op, 1))
                return;
            end
        end

        warning(['Newton''s method did not converge for ', char(string(nnz(op))), ' points (of ', char(string(numel(op))), ')']);
    end
    
    phi = w;
    w = w/C;
    d = d/C;
    
%% Large radius: Points with absolute value at least 9/8 times that of the tip
    a0 = 2*p - 1;
    a1 = p*(1 - p) / 2;
    a2 = (1 - 2*p)*p*(1 - p) / 3;
    % Iterate over z/w, instead of z
    ff = @(z, w) f(z .* w) ./ w - 1;
    dff = @(z, w) df(z .* w);
    
    lrCond = abs(w/d) > 9/8;
    phi(lrCond) = w(lrCond) .* newton(@(z, w) ff(z, w), @(z, w)(dff(z, w)), w(lrCond) + (a0 + (a1 + (a2 ./ w(lrCond))) ./ w(lrCond)), w(lrCond), false);

%% Near tip: Points within a small radius around the tip
    rota = exp(-1i*p*pi);
    k  = @(u) sqrt( rota*(u - d) );
    dk = @(u) rota / 2 ./ sqrt( rota*(u - d) );

    ntCond = ~lrCond & (abs(w - d) < 0.25*imag(d));
    phi(ntCond) = newton(@(z, w) (k(f(z)) - w), @(z, ~) (df(z) .* dk(f(z))), w(ntCond), k(w(ntCond)), true);
    
%% Right sector: Points with argument less than p*pi
    pwr = 1/p;
    fp  = @(z) (z - p) .* (z + 1 - p).^(pwr - 1);
    dfp = @(z) (z + 1 - p) .^ (pwr - 1) + (z - p) * (pwr - 1) .* (z + 1 - p) .^ (pwr - 2);
    
    rsCond = ~lrCond & ~ntCond & (angle(w) < p*pi);
    phi(rsCond) = newton(@(z, w) ( fp(z) - w ), @(z, ~) dfp(z), w(rsCond), w(rsCond).^pwr, true);
    
%% Left sector: Points with argument more than p*pi
    pwr  = 1/(1-p);
    fpp  = @(z) (z-p).^(p*pwr) .* (z+1-p);
    dfpp = @(z) p*pwr*(z-p).^(p*pwr-1) .* (z+1-p)  +  (z-p).^(p*pwr);
    
    lsCond = ~lrCond & ~ntCond & ~rsCond;
    phi(lsCond) = newton(@(z, w) ( fpp(z) - w ), @(z, ~) dfpp(z), w(lsCond), w(lsCond).^pwr, true);
    
%% Infinity maps to infinity
    phi(isinf(w)) = Inf;
    
end
