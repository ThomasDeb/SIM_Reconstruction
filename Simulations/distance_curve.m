function [d] = distance_curve(point, x_pp, y_pp)
coefs_x = x_pp.coefs; coefs_y = y_pp.coefs;
d = inf;
for i = 1 : x_pp.pieces
    temp1 = coefs_x(i, :); temp1(end) = temp1(end) - point(1);
    temp2 = coefs_y(i, :); temp2(end) = temp2(end) - point(2);
    dist = conv(temp1, temp1) + conv(temp2, temp2);
    r = roots(polyder(dist)); 
    r = r(imag(r)==0);
    % Polynomial coefs are given between 0 and x_pp.breaks(i+1) - x_pp.breaks(i)
    r(r < 0) = [];
    r(r > (x_pp.breaks(i+1) - x_pp.breaks(i))) = [];
    r = [r; 0; x_pp.breaks(i+1) - x_pp.breaks(i)];
    d = min(d, min(polyval(dist, r)));
end
d = sqrt(d);
end