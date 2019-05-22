function phi = phi_zip_fwd_np1( z, zeta1, zeta2 )
%PHI_ZIP_FWD_NP1 Forward Zipper phi_{n+1}
%   A component function of a Zipper-like algorithm.

    if (isnan(zeta1) || isnan(zeta2))
        error('Zipper error: Invalid parameter in forward 2-point final map');
    end
    
    phi = uhpaut(zeta2, z);
    d = uhpaut(zeta2, zeta1);
    
    alpha = pi - angle(d);
    phi = exp(-1i*(pi - alpha)) * phi;
    phi = phi .^ (pi / alpha);
    
end

