function [] = slit(obj)
%ZIPPIT.SLIT Set up conformal map as given by geodesic algorithm
%   Performs slit algorithm on polygon vertices to determine the parameters
%   zn, and sets up the corresponding forward and inverse for direct access.

    % Set up component functions
    % (phi_1 and phi_{n+1} are the same as geod's)
    phi_fwd_1 = @phi_geod_fwd_1;
    phi_fwd_i = @phi_slit_fwd_i;
    phi_fwd_np1 = @phi_geod_fwd_np1;

    phi_inv_1 = @phi_geod_inv_1;
    phi_inv_i = @phi_slit_inv_i;
    phi_inv_np1 = @phi_geod_inv_np1;


    % Solve parameters zn of functions phi_i
    obj.zn = [obj.polygon; obj.polygon(1)];
    obj.zn(3:end) = phi_fwd_1(obj.zn(3:end), obj.zn(1), obj.zn(2));
    for x = 3:obj.nv
        obj.zn(x+1:end) = phi_fwd_i(obj.zn(x+1:end), obj.zn(x));
    end


    % Define forward and inverse
    function w = forward(z)
        w = phi_fwd_1(z, obj.zn(1), obj.zn(2));
        for y=3:obj.nv
            w = phi_fwd_i(w, obj.zn(y));
        end
        w = phi_fwd_np1(w, obj.zn(end));
    end
    obj.forward = @forward;

    function z = inverse(w)
        z = phi_inv_np1(w, obj.zn(end));
        for y=obj.nv:-1:3
            z = phi_inv_i(z, obj.zn(y));
        end
        z = phi_inv_1(z, obj.zn(1), obj.zn(2));
    end
    obj.inverse = @inverse;

end
