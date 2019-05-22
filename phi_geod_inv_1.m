function phi = phi_geod_inv_1( z, z0, z1 )
%PHI_GEOD_INV_1 Inverse geodesic/slit phi_1
%   A component function of a Zipper-like algorithm.

    if (any(isnan([z0, z1])) || any(isinf([z0, z1])) ~= 0 || numel(z0) ~= 1 || numel(z1) ~= 1)
        error('Zipper error: Invalid parameter in inverse 2-point initial map');
    end

    sq = z.^2;
    phi = (z1 + z0*sq) ./ (sq + 1);
    

end

