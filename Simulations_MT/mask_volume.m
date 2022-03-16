function [mask] = mask_volume(sz, num_pts, sigma, radius, num_mt, num_pts_t, sigma_t, ratio_moving_pts)
%N = 256; num_t = 9; sz = [N, N, num_t];
szIm = sz(1:2); num_t = sz(3);

x_pts = zeros(num_pts, num_mt, num_t);
y_pts = zeros(num_pts, num_mt, num_t);
t_pts = round(linspace(1, num_t, num_pts_t));

% Figure
figure; num_pts_plot = 200; param_grid = linspace(0, 1, num_pts_plot);
mask = zeros(sz, 'logical');
mask_tmp = zeros(szIm, 'logical');
for i = 1 : num_mt
    [x_pts(:, i, 1), y_pts(:, i, 1), x_pp, y_pp] = generate_pts(num_pts, sigma, szIm);
    %plot(fnval(x_pp,param_grid), fnval(y_pp,param_grid)); hold on;
    %scatter(x_pts(:, i), y_pts(:, i));
    mask_tmp(generate_mask(x_pp, y_pp, radius, szIm)) = 1;
end
xlim([0 sz(1)-1]); ylim([0 sz(2)-1]);
mask(:,:,1) = mask_tmp;

for t = 2 : num_pts_t
    [x_pts(:,:,t), y_pts(:,:,t), x_pp, y_pp] = move_pts_t(x_pts(:,:,t-1), ...
        y_pts(:,:,t-1), sigma_t*(num_t-1)/(num_pts_t-1), ratio_moving_pts, szIm);
    for tt = t_pts(t-1)+1 : t_pts(t)
        mask_tmp = zeros(szIm, 'logical');
        alph = (tt - t_pts(t-1)) / (t_pts(t) - t_pts(t-1)); % Coeff for linear interpolation btw points
        for i = 1 : num_mt
            x = x_pts(:,i,t-1) +  alph * (x_pts(:,i,t) - x_pts(:,i,t-1));
            y = y_pts(:,i,t-1) +  alph * (y_pts(:,i,t) - y_pts(:,i,t-1));
            [x_pp, y_pp] = points2curve(x, y);
            if ~is_valid_curve(x_pp, y_pp, szIm)
                error('WTF');
            end
            %scatter(x_pts(:, i), y_pts(:, i));
            mask_tmp(generate_mask(x_pp, y_pp, radius, szIm)) = 1;
        end
        mask(:,:,tt) = mask_tmp;
    end
end
