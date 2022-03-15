function [mask] = mask_curve(x_pp, y_pp, sz)
mask = zeros(sz, 'logical');
grid = linspace(0, 1, numel(mask));
x = fnval(x_pp, grid); y = fnval(y_pp, grid);
if ((round(min(x))+1) < 1) || (round(max(x))+1 > sz(1)) || (round(min(y))+1 < 1) || (round(max(y))+1 > sz(2))
1;
end
mask(sub2ind(size(mask), round(x)+1, round(y)+1)) = 1;
end