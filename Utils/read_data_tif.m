function [data_vol] = read_data_tif(datapath)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

V = double(tiffreadVolume(datapath));
sz = size(V);
sz_out = [sz(1)/3, sz(2)/3, sz(3)*9];
data_vol = zeros(sz_out);
k = 1;
for slice = 1 : sz(3)
    for i = 0:2
        for j = 0:2
            data_vol(:,:,k) = V(i*sz_out(1) + (1:sz_out(1)), j*sz_out(2) + (1:sz_out(2)), slice);
            k = k + 1;
        end
    end
end
end