function phi = phi_geod_fwd_1( z, z0, z1 )
%PHI_GEOD_FWD_1 Forward geodesic/slit phi_1
%   A component function of a Zipper-like algorithm.

    if (any(isnan([z0, z1])) || any(isinf([z0, z1])) ~= 0 || numel(z0) ~= 1 || numel(z1) ~= 1)
        error('Zipper error: Invalid parameter');
    end

    phi = 1i*sqrt((z - z1) ./ (z - z0));
    
end

