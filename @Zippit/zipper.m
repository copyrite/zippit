function [] = zipper(obj)
%ZIPPIT.ZIPPER Set up conformal map as given by Zipper algorithm
%   Performs Zipper algorithm on polygon vertices to determine the parameters
%   zn, and sets up the corresponding forward and inverse for direct access.

    if (mod(obj.nv, 2) ~= 0)
        error('Zipper error: Number of polygon vertices is odd.');
    end

    % Set up component functions
    phi_fwd_1 = @phi_zip_fwd_1;
    phi_fwd_i = @phi_zip_fwd_i;
    phi_fwd_np1 = @phi_zip_fwd_np1;

    phi_inv_1 = @phi_zip_inv_1;
    phi_inv_i = @phi_zip_inv_i;
    phi_inv_np1 = @phi_zip_inv_np1;


    % Solve parameters zn of functions phi_i
    obj.zn = [obj.polygon; obj.polygon(1)];
    obj.zn(4:end) = phi_fwd_1(obj.zn(4:end), obj.zn(1), obj.zn(2), obj.zn(3));
    for x = 4:2:obj.nv-2
        obj.zn(x+2:end) = phi_fwd_i(obj.zn(x+2:end), obj.zn(x), obj.zn(x+1));
    end


    % Define forward and inverse
    function w = forward(z)
        w = phi_fwd_1(z, obj.zn(1), obj.zn(2), obj.zn(3));
        for y=4:2:obj.nv-2
            w = phi_fwd_i(w, obj.zn(y), obj.zn(y+1));
        end
        w = phi_fwd_np1(w, obj.zn(end-1), obj.zn(end));
    end
    obj.forward = @forward;

    function z = inverse(w)
        z = phi_inv_np1(w, obj.zn(end-1), obj.zn(end));
        for y=obj.nv-2:-2:4
            z = phi_inv_i(z, obj.zn(y), obj.zn(y+1));
        end
        z = phi_inv_1(z, obj.zn(1), obj.zn(2), obj.zn(3));
    end
    obj.inverse = @inverse;

end
