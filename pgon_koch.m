function koch = pgon_koch( gen )
%PGON_KOCH Koch snowflake as a complex polygon
%   Returns a generation 'gen' Koch snowflake as a complex vector of
%   vertices.

    function new = iter(old)
        x = old;
        y = circshift(old, -1);
        new = reshape([3*x, 2*x + y, 2*x+y + exp(-1i*pi/3)*(y - x), x + 2*y].'/3, [4*numel(old), 1]);
    end

if (numel(gen) ~= 1)
    error('Erroneous parameter given to Koch snowflake');
end

koch = zeros(3*4^gen, 1);
koch(1:3) = 1i*exp((0:2) * 2i * pi/3);

for i=1:gen
    koch(1:(3*4^i)) = iter(koch(1:(3*4^(i-1))));
end

end
