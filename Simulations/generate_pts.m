function [x_n, y_n, x_pp, y_pp] = generate_pts(num_pts, sigma, sz)
x = rand([2, 1]) * sz(1);
y = rand([2, 1]) * sz(2);

%Put start and end points on border
r = randperm(4, 2);
switch r(1)
    case 1
        x(1) = 0;
    case 2
        x(1) = sz(1)-1;
    case 3
        y(1) = 0;
    case 4
        y(1) = sz(2)-1;
end
switch r(2)
    case 1
        x(end) = 0;
    case 2
        x(end) = sz(1)-1;
    case 3
        y(end) = 0;
    case 4
        y(end) = sz(2)-1;
end

v_orth = [-(y(end) - y(1)), x(end) - x(1)]; v_orth = v_orth / norm(v_orth);
x = linspace(x(1), x(end), num_pts);
y = linspace(y(1), y(end), num_pts);
sigma_temp = sigma * norm([x(1) - x(end), y(1) - y(end)]) / max(sz); % Normalize sigma depending on length
isValid = false;
while ~isValid
    x_n = x; y_n = y;
    for i = 2 : num_pts-1
        r = randn;
        x_n(i) = x(i) + sigma_temp * r * v_orth(1);
        y_n(i) = y(i) + sigma_temp * r * v_orth(2);
        while (x_n(i) < 0 || x_n(i) > (sz(1)-1) || y_n(i) < 0 || y_n(i) > (sz(2)-1))
            r = randn;
            x_n(i) = x(i) + sigma_temp * r * v_orth(1);
            y_n(i) = y(i) + sigma_temp * r * v_orth(2);
        end
        %     y_n(i) = y(i) + sigma * randn;
        %     while (y_n(i) < 0 || y_n(i) > (sz(2)-1))
        %         y_n(i) = y(i) + sigma * randn;
        %     end
        [x_pp, y_pp] = points2curve(x_n, y_n);
        isValid = is_valid_curve(x_pp, y_pp, sz);
    end
end