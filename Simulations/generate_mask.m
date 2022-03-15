function [mask] = generate_mask(x_pp, y_pp, radius, sz)
mask = zeros(sz, 'logical');
D = bwdist(mask_curve(x_pp, y_pp, sz));
[ind_x, ind_y] = find(D <= (radius+2));
for j = 1 : length(ind_x)
    d = distance_curve([ind_x(j), ind_y(j)], x_pp, y_pp);
    mask(ind_x(j), ind_y(j)) = (d <= radius);
end
end