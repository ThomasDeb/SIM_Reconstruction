function [x_pp, y_pp] = points2curve(x, y)
num_pts = length(x);
param_knots = linspace(0, 1, num_pts);
x_pp = csapi(param_knots, x);
y_pp = csapi(param_knots, y);
end