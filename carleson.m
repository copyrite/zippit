function grid = carleson( gen )
%CARLESON Vertices of a Carleson grid
%   Returns the vertices of a generation 'gen' Carleson grid. The grid has
%   2^gen squares on the outermost rim. gen must be at least 2.

if (numel(gen) ~= 1 || (gen < 2))
    error('Erroneous parameter given to Carleson grid');
end

grid = zeros(2^(gen+1) - 4, 1);
nGrid = 0;
for x=2:gen
    grid(nGrid + (1:2^x)) = (1 - (pi/2)^-x)*exp(2i*pi*(1:2^x)' / 2^x);
    nGrid = nGrid + 2^x;
end

end

