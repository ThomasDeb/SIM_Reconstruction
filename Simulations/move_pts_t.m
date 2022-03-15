function [x_pts_new, y_pts_new, x_pp, y_pp] = move_pts_t(x_pts, y_pts, sigma, ratio_moving_pts, sz)

num_pts = numel(x_pts);
num_moving_pts = round(ratio_moving_pts * num_pts);
moving_pts = randperm(num_pts, num_moving_pts);
x_pts_new = x_pts;
y_pts_new = y_pts;

for i = 1 : num_moving_pts
    [idx_pts, idx_mt] = ind2sub(size(x_pts), moving_pts(i));
    isValid = false;
    while ~isValid
        x_pts_new(idx_pts, idx_mt) = x_pts(idx_pts, idx_mt) + sigma * randn;
        y_pts_new(idx_pts, idx_mt) = y_pts(idx_pts, idx_mt) + sigma * randn;
        [x_pp, y_pp] = points2curve(x_pts_new(:, idx_mt), y_pts_new(:, idx_mt));
        isValid = is_valid_curve(x_pp, y_pp, sz);
    end
end

end