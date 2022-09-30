function [mask] = mask_sim(sz,orr,ns,lamb,res,bet,fc)
[X,Y]= meshgrid(1:sz(1),1:sz(2));
mask = zeros(sz(1:2));
for ii=1:length(orr)
    k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
    f_idx = k/pi*res.*sz(1:2);
    center_plus = (sz(1:2)+1)/2-f_idx;
    center_min = (sz(1:2)+1)/2+f_idx;
    mask_plus = ((X-center_plus(1)).^2 + (Y-center_plus(2)).^2) <= (fc*sz(1))^2;
    mask_min = ((X-center_min(1)).^2 + (Y-center_min(2)).^2) <= (fc*sz(1))^2;
    mask(mask_plus) = 1; mask(mask_min) = 1;
end
end