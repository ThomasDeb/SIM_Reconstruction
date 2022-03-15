N = 256; num_t = 9; sz = [N, N, num_t];
% Microtubule parameters
rerun = 0;
% Spatial parameters
num_pts = 5; sigma = 10; radius = 1; num_mt = 10;
% Temporal parameters
num_pts_t = 5; sigma_t = 1; ratio_moving_pts = 0.8;

if rerun
    im = mask_volume(sz, num_pts, sigma, radius, num_mt, num_pts_t, sigma_t, ratio_moving_pts);
    saveastiff(double(im),[expFolder,'objectUFS_MT.tif']);
    %save([expFolder,'objectUFS_MT.mat'], im)
else
    im=double(loadtiff([expFolder,'objectUFS_MT.tif']));
    %load([expFolder,'objectUFS_MT.mat'], im);
end

%figure; sliceViewer(im, 'Colormap', parula(256)); colorbar;

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
patt=zeros(sz);
[X,Y]=meshgrid(0:sz(2)-1,0:sz(1)-1);X=X*res;Y=Y*res;
i=1;
for ii=1:length(orr)
    k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
    for jj=1:length(ph)
        patt(:,:,i)=1+ a*cos(2*(k(1)*X+k(2)*Y + ph(jj)));
        i=i+1;
    end
end
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

%% Low-passed images
%lamb=488; res=32; Na=1.4; fc=2*Na/lamb*res; nl=1.518; ns=1.333; orr=[0 pi/3 2*pi/3]; bet=asin(Na/nl);     
mask_lp = mask_lowpass(sz, fc);
mask_SIM = mask_sim(sz, orr, ns, lamb, res, bet, fc);
mask_FT = zeros(sz);
im_lowpass = zeros(sz);
im_sim = zeros(size(im));
for t = 1 : size(im, 3)
    mask_FT(:,:,t) = fft2(squeeze(im(:,:,t)));
    im_lowpass(:,:,t) = ifft2(ifftshift(fftshift(fft2(im(:,:,t))) .* mask_lp));
    im_sim(:,:,t) = ifft2(ifftshift(fftshift(fft2(im(:,:,t))) .* mask_SIM));
end
%figure; sliceViewer(log(1+abs(fftshift(mask_FT))), 'Colormap', parula); colorbar; title('Ground truth FFT')
figure; sliceViewer(abs(im_lowpass), 'Colormap', parula); colorbar; title('OTF low-passed ground truth')
figure; sliceViewer(abs(im_sim), 'Colormap', parula); colorbar; title('SIM low-passed ground truth')
