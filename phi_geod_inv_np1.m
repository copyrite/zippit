function phi = phi_geod_inv_np1( z, zeta )
%PHI_GEOD_INV_NP1 Inverse geodesic/slit phi_{n+1}
%   A component function of a Zipper-like algorithm.

    if (isnan(zeta))
        error('Zipper error: Invalid parameter in inverse 1-point final map');
    end

    phi = 1i*sqrt(z);
    phi = uhpaut(-zeta, phi);
    
end

