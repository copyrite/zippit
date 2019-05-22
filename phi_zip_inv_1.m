function phi = phi_zip_inv_1( z, z0, z1, z2 )
%PHI_ZIP_INV_1 Inverse Zipper phi_1
%   A component function of a Zipper-like algorithm.

    if (any(isnan([z0, z1, z2])) || any(isinf([z0, z1, z2])) ~= 0 || numel(z0) ~= 1 || numel(z1) ~= 1 || numel(z2) ~= 1)
        error('Zipper error: Invalid parameter in inverse 3-point initial map');
    end
    
% Inverting phi_zip_fwd_1:
%   w   = sqrt( z - z2  *  z1 - z0  /  z - z0  /  z1 - z2 )
%   w^2 =      (z - z2) * (z1 - z0) / (z - z0) / (z1 - z2)
%   w^2 (z - z0) (z1 - z2)          = (z - z2) * (z1 - z0)
%   z (w^2 (z1 - z2) - (z1 - z0))   = z0 w^2 (z1 - z2) - z2 (z1 - z0)
%   z = z0 w^2 (z1 - z2) - z2 (z1 - z0) / (w^2 (z1 - z2) - (z1 - z0))

    sq = z.^2;
    phi = ( z0*sq*(z1-z2) - z2*(z1-z0) ) ./ ( sq*(z1-z2) - (z1-z0) );
    
end

