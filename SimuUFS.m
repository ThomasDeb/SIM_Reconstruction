%----------------------------------------------------
% Script to simulate SIM data. Illumination patterns are generated using:
%
%  W(x,y) = 1 + a cos(2*(kx x + ky y + ph))
%
% where a is the amplitude of the modulation, ph the phase of the
% illumination pattern and k=[kx,ky] is the wave vector given by
%   k=(2 pi ns)/lamb [cos(th), sin(th)]*sin(bet)
% with th the lateral orientation of the illumination, nl the refractive
% index of the objective medium and bet the angle between side beams and
% the optic axis (of the two interfering beams used to generate the
% pattern) which should be lower than asin(NA/nl).
%
% If the given image is 3D, the script perform a 3D acquisition and returns
% only the focal plane.
%
% Before to run it set the following papareters:
% -- General
% sav=             -> to save results
% expFolder=       -> experiment folder
% gtpath=          -> file name ground truth
%
% -- PSF
% lamb=            -> Illumination wavelength
% res=             -> Resolution
% Na=              -> Objective numerica aperture
% nl=              -> Refractive index of the objective medium (glass/oil)
% ns=              -> Refractive index of the sample medium (water)
%
% -- Patterns
% orr=             -> Patterns orientations (vector)
% ph=              -> Patterns lateral phases (vector)
% a=               -> Amplitude coefficient
% bet=             -> Angle between side beams and the optic axis (e.g. bet asin(Na/nl))
%
% -- Acquisition
% downFact=        -> Downsmpling factor (e.g. [2 2])
% photBud=         -> Photon Budget (average per voxel of ground truth)
%
% Copyright (2018) Emmanuel Soubies, emmanuel.soubies@epfl.ch
%----------------------------------------------------
addpath ./Utils
javaaddpath ./Utils/PSFGenerator.jar

%% Reading data
fprintf('Reading data .............');
im2D=double(read(Tiff(gtpath, 'r')));
fprintf(' done \n');
im2D = imresize(im2D, gtsize);
im2D = (im2D - min(im2D(:)))/ (max(im2D(:)) - min(im2D(:)));
nbPatt = length(orr)*length(ph)*num_cycles;
sz = [size(im2D), nbPatt];
bckgrnd = im2D(5,5); % Background value
im2D(1:5, :) = bckgrnd; im2D(end-4:end, :) = bckgrnd;
im2D(:, 1:5) = bckgrnd; im2D(:, end-4:end) = bckgrnd;
% Moving image
% max_displacement_x = 5; % Max displacement in total
% displacement_x = round(linspace(0, max_displacement_x, ceil((nbPatt+1)/2)));
% displacement_x = [displacement_x, flip(displacement_x(2-mod(nbPatt, 2):end-1))];
% max_displacement_y = max_displacement_x; %
% displacement_y = round(linspace(0, max_displacement_y, ceil((nbPatt+1)/2)));
% displacement_y = [displacement_y, flip(displacement_y(2-mod(nbPatt, 2):end-1))];
max_displacement_x = nbPatt; % Max displacement in total
displacement_x = round(linspace(0, max_displacement_x, nbPatt));
max_displacement_y = max_displacement_x; %
displacement_y = round(linspace(0, max_displacement_y, nbPatt));
im = ones(sz) * bckgrnd; 

for n = 1 : nbPatt
    if displacement_x(n) > 0
        intx = (1+displacement_x(n)) : sz(1);
    else
        intx = 1 : (sz(1)-displacement_x(n));
    end
    if displacement_y(n) > 0
        inty = (1+displacement_y(n)) : sz(2);
    else
        inty = 1 : (sz(2)-displacement_y(n));
    end
    im(intx,inty,n) = im2D(intx-abs(displacement_x(n)), inty-abs(displacement_y(n)));
end
saveastiff(double(im),[expFolder,'objectUFS.tif']);


%% PSF Generation
fprintf('PSF Generation ...........');
fc=2*Na/lamb*res;
ll=linspace(-0.5,0,sz(1)/2+1);
lr=linspace(0,0.5,sz(1)/2);
[X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
[th,rho]=cart2pol(X,Y);
OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
figure;subplot(1,2,1);imagesc((fftshift(OTF))); axis image; title('OTF');
colormap(fire(200));viscircles(floor(sz(1:2)/2)+1,fc*sz(1));
psf=real(fftshift(ifft2(OTF)));
OTF = fft2(fftshift(psf)); % Bug fix ?
subplot(1,2,2);imagesc(psf); axis image; title('PSF'); caxis([0 0.01*max(psf(:))]);
fprintf(' done \n');

%% Patterns Generation
fprintf('Patterns Generation ......');
patt=zeros([sz(1:2), length(orr)*length(ph)]);
[X,Y]=meshgrid(0:sz(2)-1,0:sz(1)-1);X=X*res;Y=Y*res;
i=1;
for ii=1:length(orr)
    k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
    for jj=1:length(ph)
        patt(:,:,i)=1+ a*cos(2*(k(1)*X+k(2)*Y + ph(jj)));
        i=i+1;
    end
end
patt = repmat(patt, [1, 1, num_cycles]);
nbPatt=size(patt,3); % Normalization such that the mean of each pattern is 1/#Patterns
for ii=1:nbPatt
    tmp=patt(:,:,ii);
    patt(:,:,ii)=patt(:,:,ii)/(mean(tmp(:))*nbPatt);
end
nn=sum(patt.^2,3);                    
patt=patt/(sqrt(max(nn(:)))); % Normalize patterns for simplicity (see Manu's paper)
figure;subplot(1,2,1);imagesc(patt(:,:,1)); axis image; title('Example pattern');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(patt(:,:,1)))))); axis image; title('Example pattern FFT');
viscircles(floor(sz(1:2)/2)+1,fc*sz(1));
fprintf(' done \n');

%% Image Noiseless Acquisition
fprintf('Acquisition simulation ...');
% - LinOp Downsampling and integration over camera pixels

S=LinOpDownsample(sz(1:2),downFact);
% - LinOpConv (PSF)
%OTF=Sfft(fftshift(fftshift(psf(:,:),1),2),3);
H=LinOpConv(OTF,1,[1 2]);
% - Acquisition
acqNoNoise=zeros([S.sizeout,size(patt,3)]);
fprintf(' Pattern # ');
for i=1:size(patt,3)
    fprintf([num2str(i),'..']);
    Si = LinOpSelectorPlanes(sz,3,i,true);
    Di = LinOpDiag(sz(1:2),patt(:,:,i));
    acqNoNoise(:,:,i)=S*H*Di*Si*im;
end
fprintf(' done \n');

%% Add noise and Save
for ii=1:length(photBud)
    % - Add Noise
    acq=acqNoNoise;
    if photBud>0
        tmp=sum(acqNoNoise,3);
        factor = photBud(ii)./mean(tmp(:)) ;
        acqNoNoise = acqNoNoise.* factor;
        acq = random('Poisson',acqNoNoise);
        im = im.*factor;
    end
    SNR=20*log10(norm(acqNoNoise(:))/norm(acq(:)-acqNoNoise(:)));
    disp(['SNR = ',num2str(SNR),' dB']);
    
    % - Save
    if sav
        saveastiff(double(acqNoNoise),[expFolder,'AcqDataNoiseless.tif']);
        saveastiff(double(log(1+abs(fftshift(fftshift(Sfft(acqNoNoise,3),1),2)))),[expFolder,'AcqDataNoiseless-FFT.tif']);
        saveastiff(double(acq),[expFolder,'AcqData.tif']);
        saveastiff(double(log(1+abs(fftshift(fftshift(Sfft(acq,3),1),2)))),[expFolder,'AcqData-FFT.tif']);
        saveastiff(double(sum(acq,3)),[expFolder,'WFData.tif']);
        saveastiff(double(log(1+abs(fftshift(fft2(sum(acq,3)))))),[expFolder,'WFData-FFT.tif']);
        saveastiff(double(sum(acqNoNoise,3)),[expFolder,'WFDataNoiseless.tif']);
        saveastiff(double(log(1+abs(fftshift(fft2(sum(acqNoNoise,3)))))),[expFolder,'WFDataNoiseless-FFT.tif']);
    end
end
save([expFolder,'PSF'],'psf');
saveastiff(double(patt),[expFolder,'patterns.tif']);

figure;subplot(1,2,1);imagesc(acq(:,:,1)); axis image; title('Example Acquired image');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(acq(:,:,1)))))); axis image; title('Example Acquired image FFT');
figure;subplot(1,2,1);imagesc(sum(acq,3)); axis image; title('Widefield image');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(sum(acq,3)))))); axis image; title('Widefield image FFT');

%% Show all sampled frequencies and low-passed image

figure; imagesc(log(1+abs(fftshift(fftn(im2D))))); axis image; title('Object FFT');
[X,Y]=meshgrid(1:sz(1),1:sz(2));
mask_lowpass = ((X-(sz(1)+1)/2).^2 + (Y-(sz(2)+1)/2).^2) <= (fc*sz(1))^2;
mask_sim = zeros(sz(1:2));
for ii=1:length(orr)
    k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
    f_idx = k/pi*res.*sz(1:2);
    center_plus = (sz(1:2)+1)/2-f_idx;
    center_min = (sz(1:2)+1)/2+f_idx;
    viscircles(center_plus,fc*sz(1), 'color', 'b');
    viscircles(center_min,fc*sz(1), 'color', 'b');
    mask_plus = ((X-center_plus(1)).^2 + (Y-center_plus(2)).^2) <= (fc*sz(1))^2;
    mask_min = ((X-center_min(1)).^2 + (Y-center_min(2)).^2) <= (fc*sz(1))^2;
    mask_sim(mask_plus) = 1; mask_sim(mask_min) = 1;
end
viscircles(floor(sz(1:2)/2)+1,fc*sz(1));

im2D_lowpass = ifft2(ifftshift(fftshift(fft2(im2D)) .* mask_lowpass));
im2D_sim = ifft2(ifftshift(fftshift(fft2(im2D)) .* mask_sim));
figure;imagesc(abs(im2D_lowpass)); axis image; title('Lowpass object');
figure;imagesc(abs(im2D_sim)); axis image; title('Lowpass object (9 circles)');



