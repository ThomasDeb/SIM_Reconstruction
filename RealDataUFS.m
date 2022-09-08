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

%% Load data
fprintf('Reading data .............');
y=double(tiffreadVolume([dataFolder, dataFile, '/', dataFile, '_stack.tif'])); y=y/max(y(:));
szacq=size(y); sz = [szacq(1:2)*2, szacq(3)];
fprintf(' done \n');

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

[patt,orr,phOff,k0]=PatternsFromFairSimTxt([dataFolder, dataFile, '/patterns_a=01.txt']);
nbPatt = size(patt,3);

nn=sum(patt.^2,3);                    
patt=patt/(sqrt(max(nn(:)))); % Normalize patterns for simplicity (see Manu's paper)

figure;subplot(1,2,1);imagesc(patt(:,:,1)); axis image; title('Example pattern');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(patt(:,:,1)))))); axis image; title('Example pattern FFT');
viscircles(floor(sz(1:2)/2)+1,fc*sz(1));
fprintf(' done \n');

%% Patterns on coarse grid
    a = 0.1;
    patt_coarse=zeros([szacq(1:2),nbPatt]);
    phShift=linspace(0,2*pi,length(phOff)+1); phShift=phShift(1:end-1); 
    [X,Y]=meshgrid(0:szacq(1)-1,0:szacq(2)-1);
    it=1;
    for ii=1:length(orr)
        k=2*pi*[cos(orr(ii)), sin(orr(ii))]*(k0(ii))/(szacq(1)*2);
        for jj=1:length(phOff)
            patt_coarse(:,:,it) = 1 + a*cos(2*(k(1)*X+k(2)*Y + phOff(ii)/2+ phShift(jj)/2));
            it=it+1;
        end
    end
    
    y_no_patt = zeros(size(y));
    for t = 1 : size(y, 3)
        y_no_patt(:,:,t) = y(:,:,t) ./ patt_coarse(:,:,1+mod(t-1, size(patt, 3)));
    end
    figure; sliceViewer(y_no_patt, 'Colormap', parula); title('Data without illumination patterns'); colorbar;


%% Save
save([expFolder,'PSF'],'psf');
saveastiff(double(patt),[expFolder,'patterns.tif']);


