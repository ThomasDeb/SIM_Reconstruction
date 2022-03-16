function [mask] = mask_lowpass(sz, fc)
[X,Y]= meshgrid(1:sz(1),1:sz(2));
mask = ((X-(sz(1)+1)/2).^2 + (Y-(sz(2)+1)/2).^2) <= (fc*sz(1))^2;
end