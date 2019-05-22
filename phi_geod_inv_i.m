function phi = phi_geod_inv_i( z, a )
%PHI_GEOD_INV_I Inverse geodesic phi_i
%   A component function of a Zipper-like algorithm.

    if (numel(a) ~= 1 || isnan(a) || isinf(a))
        error('Zipper error: Invalid parameter');
    end
    
    b = abs(a)^2 / real(a);
    c = abs(a)^2 / imag(a);

    phi = z .* sqrt(1 - c^2 * (z.^-2));
    
    phi = uhpaut(-b, phi);
    
end

