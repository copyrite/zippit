function phi = phi_geod_fwd_i( z, a )
%PHI_GEOD_FWD_I Forward geodesic phi_i
%   A component function of a Zipper-like algorithm.

    if (numel(a) ~= 1 || isnan(a) || isinf(a))
        error('Zipper error: Invalid parameter');
    end

    b = abs(a)^2 / real(a);
    c = abs(a)^2 / imag(a);
    
    phi = uhpaut(b, z);
    phi = phi .* sqrt(1 + c^2 * (phi.^-2));
    
end

