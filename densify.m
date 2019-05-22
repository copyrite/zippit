function p = densify(pgon, num)
%DENSIFY Add vertices to a polygon boundary
%
%   DENSIFY (pgon, num) returns a polygon, represented as a vector of
%   complex numbers, that has all the vertices of 'pgon' but each line
%   has 'num' amount of uniformly added new vertices.
%
%   Every (num+1)th vertex of the vector is from the original polygon,
%   starting from index 1.

intv = (0:num)/(num+1);
p = reshape((pgon + ([pgon(2:end); pgon(1)] - pgon)*intv).', [numel(pgon) * (num+1), 1]);

end