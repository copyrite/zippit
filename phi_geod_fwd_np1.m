function phi = phi_geod_fwd_np1( z, zeta )
%PHI_GEOD_FWD_NP1 Forward geodesic/slit phi_{n+1}
%   A component function of a Zipper-like algorithm.

    if (isnan(zeta))
        error('Zipper error: Invalid parameter in forward 1-point final map');
    end

    phi = -(uhpaut(zeta, z).^2);

end

