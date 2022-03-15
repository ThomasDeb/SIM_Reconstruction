function [b] = is_valid_curve(x_pp, y_pp, sz)
xmin = fnmin(x_pp, [0, 1]);
ymin = fnmin(y_pp, [0, 1]);
x_pp.coefs = - x_pp.coefs;
y_pp.coefs = - y_pp.coefs;
xmax = - fnmin(x_pp, [0, 1]);
ymax = - fnmin(y_pp, [0, 1]);

b = (xmin >= 0) && (ymin >= 0) && (xmax <= (sz(1)-1)) && (ymax <= (sz(2)-1));
end