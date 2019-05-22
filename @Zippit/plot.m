function [fig, cgrid] = plot( zip, ctr )
% Zippit.PLOT Plots a Zipper-like map

fig = figure;

cgrid = zip.inverse(udisk2uhp(carleson(12), zip.forward(ctr)));

cla;
plot(cgrid, '.')
axis equal;

end
