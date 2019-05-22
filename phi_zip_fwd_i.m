function phi = phi_zip_fwd_i( z, z1, z2 )
%PHI_ZIP_FWD_I Forward Zipper phi_i
%   A component function of a Zipper-like algorithm.

    if (numel(z1) ~= 1 || numel(z2) ~= 1 || any(isnan([z1, z2])) || any(isinf([z1, z2])))
        error('Zipper error: Invalid parameter in forward Zipper map');
    end
    
    % Zipper's phi_i is composed of an upper half-plane automorphism and a
    % slit map. Find the parameter of the former, then apply both.
    ctr = ([1 1i] * (2 * [real(z1) imag(z1); real(z2) imag(z2)] \ [abs(z1)^2; abs(z2)^2]));
    b = 2 * real(ctr);
    
    phi = uhpaut(b, z);
    d = uhpaut(b, z2);
    
    phi = phi_slit_fwd_i(phi, d);
    
end

