im_name = 'test_SIM026';
t_old = Tiff([im_name, '.tif'],'r');
vol9im = tiffreadVolume([im_name, '.tif']);
%imageData = read(im);
sz = size(vol9im);
figure; sliceViewer(vol9im, 'Colormap', parula);
imsz = [sz(1)/3, sz(2)/3];
vol = zeros([imsz, sz(3)*9]);
t_idx = 1;
for t = 1 : size(vol9im, 3)
    for i = 1 : 3
        for j = 1 : 3
            vol(:,:,t_idx) = vol9im(((i-1)*imsz(1)+1):i*imsz(1), ((j-1)*imsz(2)+1):j*imsz(2), t);
            t_idx = t_idx + 1;
        end
    end
end
figure; sliceViewer(vol, 'Colormap', parula);

%% View FFTs
fftvol = zeros(size(vol));
for i = 1 : size(vol, 3)
    fftvol(:,:,i) = log(1+abs(fftshift(fft2(vol(:,:,i)))));
end
figure; sliceViewer(fftvol, 'Colormap', parula); colorbar;

%% Save Tiff file
voluint16 = uint16(vol);
t = Tiff([im_name, '_stack.tif'],'w');
tagstruct.ImageLength = size(voluint16,1);
tagstruct.ImageWidth = size(voluint16,2);
tagstruct.SamplesPerPixel = getTag(t_old,'SamplesPerPixel');
tagstruct.Photometric = getTag(t_old,'Photometric');
tagstruct.PlanarConfiguration = getTag(t_old,'PlanarConfiguration');
tagstruct.BitsPerSample = getTag(t_old,'BitsPerSample');
setTag(t,tagstruct);
write(t, squeeze(voluint16(:,:,1)));
for i = 2 : size(voluint16, 3)
    writeDirectory(t);
    setTag(t,tagstruct);
    write(t, squeeze(voluint16(:,:,i)));
end
close(t);

%% Save Tiff file by time frame
voluint16 = uint16(vol);
tagstruct.ImageLength = size(voluint16,1);
tagstruct.ImageWidth = size(voluint16,2);
tagstruct.SamplesPerPixel = getTag(t_old,'SamplesPerPixel');
tagstruct.Photometric = getTag(t_old,'Photometric');
tagstruct.PlanarConfiguration = getTag(t_old,'PlanarConfiguration');
tagstruct.BitsPerSample = getTag(t_old,'BitsPerSample');

mkdir(im_name);
for t = 1 : size(vol9im, 3)
    ti = Tiff(['./', im_name, '/', num2str(t), '.tif'],'w');
    setTag(ti,tagstruct);
    write(ti, squeeze(voluint16(:,:,9*(t-1)+1)));
    for i = 2 : 9
        writeDirectory(ti);
        setTag(ti,tagstruct);
        write(ti, squeeze(voluint16(:,:,9*(t-1)+i)));
    end
    close(ti);
end

