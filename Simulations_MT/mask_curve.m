function [mask] = mask_curve(x_pp, y_pp, sz)
mask = zeros(sz, 'logical');
grid = linspace(0, 1, prod(sz));
x = fnval(x_pp, grid); y = fnval(y_pp, grid);
mask(sub2ind(size(mask), round(x)+1, round(y)+1)) = 1;
end