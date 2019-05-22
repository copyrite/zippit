function phi = phi_zip_fwd_1( z, z0, z1, z2 )
%PHI_ZIP_FWD_1 Forward Zipper phi_1
%   A component function of a Zipper-like algorithm.

    if (any([numel(z0), numel(z1), numel(z2)] ~= 1) || any(isnan([z0, z1, z2])) || any(isinf([z0, z1, z2])) ~= 0)
        error('Zipper error: Invalid parameter in forward 3-point initial map');
    end

    phi = 1i*sqrt( (z1 - z0) ./ (z1 - z2) .* (z2 - z) ./ (z - z0) );
    
end

