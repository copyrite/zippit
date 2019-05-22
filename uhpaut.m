function u = uhpaut(b, z)
%UHPAUT Upper half-plane automorphism
%   Automorphism of the upper half-plane, parameterized by the real number b.
%   Conformally maps the upper half-plane onto itself by shifting the real axis
%   as follows:
%
%   0 maps to b
%   b maps to Inf
%   Inf maps to -b
%
%   The Möbius transformation having these properties is z -> z / (1 - z/b).
%
%   In addition to reals, b is allowed to be Inf. In this case UHPAUT is an
%   identity map on z.

    if (numel(b) ~= 1 || isnan(b))
        error('Upper half-plane automorphism error: Invalid parameter')
    end
    
    u = z;
    if (~isinf(b))
        u = u ./ (1 - u/b);
    end
    
    u(isinf(z)) = -b;

end