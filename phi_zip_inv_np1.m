function phi = phi_zip_inv_np1( z, zeta1, zeta2 )
%PHI_ZIP_INV_NP1 Inverse Zipper phi_{n+1}
%   A component function of a Zipper-like algorithm.

    if (isnan(zeta1))
        error('Zipper error: Invalid parameter in inverse 2-point final map');
    end
    
    d = uhpaut(zeta2, zeta1);
    
    alpha = pi - angle(d);
    phi = z .^ (alpha/pi);
    phi = exp(1i*(pi - alpha)) * phi;
    
    phi = uhpaut(-zeta2, phi);
    
end

