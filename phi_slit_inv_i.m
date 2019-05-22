function phi = phi_slit_inv_i( z, a )
%PHI_SLIT_INV_I Inverse slit phi_i
%   A component function of a Zipper-like algorithm.

        if (numel(a) ~= 1 || isnan(a) || isinf(a))
            error('Zipper error: Invalid parameter in inverse slit map');
        end
        
        p = angle(a) / pi;
        C = abs(a) / (p^p * (1-p)^(1-p));
        
        phi = C * (z - p).^p .* (z + 1 - p).^(1 - p);
end
