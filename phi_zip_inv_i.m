function phi = phi_zip_inv_i( z, z1, z2 )
%PHI_ZIP_INV_I Inverse Zipper phi_i
%   A component function of a Zipper-like algorithm.

    if (numel(z1) ~= 1 || numel(z2) ~= 1 || any(isnan([z1, z2])) || any(isinf([z1, z2])))
        error('Zipper error: Invalid parameter in inverse Zipper map');
    end
    
    % Find the necessary parameters b and d first.
    ctr = ([1 1i] * (2 * [real(z1) imag(z1); real(z2) imag(z2)] \ [abs(z1)^2; abs(z2)^2]));
    b = 2 * real(ctr);
    d = uhpaut(b, z2);

    phi = phi_slit_inv_i(z, d);
    phi = uhpaut(-b, phi);
    
end

