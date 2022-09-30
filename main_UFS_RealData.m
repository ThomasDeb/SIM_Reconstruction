clear; close all;
setGlobalBioImPath;
addpath('./Utils/');
rng(1);

%% Generate PSF and patterns
% -- Paths
expFolder='./RealData/';   % experiment folder
dataFolder = './Data_Arne/goodPixelSize/Nice/';
dataFile = 'test_SIM026';
sav=1;                           % to save results

% -- PSF
lamb= 605;          % Illumination wavelength (EmWavelength)
res= 32;            % Resolution (Voxel size / 2)
Na= 1.49;           % Objective numerical aperture (Numerical Aperture)
nl= 1.515;          % Refractive index of the objective medium (glass/oil) (Refractive Index)
ns= 1.333;          % Refractive index of the sample medium (water) (????)

% -- Acquisition
downFact = [2 2];

% -- Run Simulator
disp('######## Generate SIM data ########');
run './RealDataUFS.m'
disp('###################################');

%% Reconstruction
% -- Paths
basedir='./RealData/';
psfname=[basedir,'PSF'];               % File name psf
pattname=[basedir,'patterns.tif'];     % File name patterns
gtname=[];           % File name ground truth (if no let empty)
outFolder=[basedir];                    % Folder to save results
sav = true;                                  % Boolean if true save the result
do_noTV = false;
% -- Data
valback=0;       % Background value

%% Reconstruction
disp('########## Reconstruction #########');

global isGPU

%% Reading data
% -- PSF
if strcmp(psfname(end-2:end),'tif')
    psf=double(loadtiff(psfname));%psf=psf/sum(psf(:));
else
    load(psfname);
end
% -- Patterns
if ~isempty(pattname)
    patt=double(loadtiff(pattname));
    if min(patt(:)) < 0, patt=patt-min(patt(:)); end                       % Ensure that patterns are positive]
    nbPatt=size(patt,3);
end
% -- Data
nbTime = 2;
y = y(:,:,1:nbPatt*nbTime)-valback; %maxy=max(y(:));y=y/maxy;              % Load and normalize data
y = imresize3(y,[512 512 size(y, 3)]);
figure;sliceViewer(y, 'Colormap', parula); title('Data'); colorbar; caxis([min(y(:)), max(y(:))]);
% -- Ground Truth
if ~isempty(gtname)                                                        % Load ground truth if available
    gt=double(loadtiff(gtname))-valback;
else
    gt=[];
end
% -- Some parameters
szUp=size(psf);
sz=[szUp, size(y, 3)]; %input size

% -- Objective Function
lamb_TV=3*1e-3 * mean(y(:));       % Hyperparameter (can be an array to loop)
lamb_hess=lamb_TV*1;       % Hyperparameter (can be an array to loop)
symPsf=1;        % Boolean true if psf is symmetric

% -- SIM reconstruction
maxIt = 100;               % Max iterations
ItUpOut = round(maxIt/10);  % Iterations between to call to OutputOpti
rhoDTNN = 1e-4;            % rho parameter (ADMM) associated to data term
rho_TV = lamb_TV*0.1;             % rho parameter (ADMM) associated to TV term (must be greater or equal than rhoDTNN if iterCG=0)
rho_hess = lamb_hess*0.01;             % rho parameter (ADMM) associated to Hessian-Schatten term (must be greater or equal than rhoDTNN if iterCG=0)
valId = 2;                  % Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 3)
periodic_TV = false;

%% Conversion CPU/GPU is necessary
psf=gpuCpuConverter(psf);
y=gpuCpuConverter(y);

%% Common operators and costs
% -- Data term
%psf=psf/sum(psf(:));
otf=fft2(fftshift(psf));
H=LinOpConv(otf, 1, 1:2);
S=LinOpDownsample(szUp,downFact);
Fn={};Hn={};F0=[];
if ~isempty(pattname)
    nn=sum(patt.^2,3);                      % Trick to perform the last ADMM step in Fourier
    patt=patt/(sqrt(max(nn(:))));
    for t = 1 : size(y, 3)
        St = LinOpSelectorPlanes(sz,3,t,true);
        L2=CostL2([],y(:,:,t));             % L2 cost function
        Fn{t}=(L2*(S*H));
        Fn{t}.doPrecomputation=true;
        Hn{t}=LinOpDiag(H.sizein,patt(:,:,1+mod(t-1, size(patt, 3))))*St;
    end
else
    L2=CostL2([],y);             % L2 cost function
    LS=L2*(S*H);
    LS.doPrecomputation=true;
end
% -- Regularization
% TV regularization in time
if periodic_TV
    Op_TV = LinOpGrad(sz, 3, 'circular');      % TV regularizer: Operator Gradient
    F_TV = CostL1(Op_TV.sizeout);                % TV regularizer: L1 norm (gradient of size 1)
else
    Op_TV = LinOpIdentity(sz);
    F_TV = CostTV1D(sz, 3);
end
Op_hess=LinOpHess(sz,'circular',[1 2]);  % Hessian-Shatten: Hessian Operator
F_hess=CostMixNormSchatt1(Op_hess.sizeout,1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
% -- Non-Negativity constraint
pos=CostNonNeg(sz);
assert(valId >= 1,'valId must be greater than 1');
OpD=LinOpDiag(sz,sqrt(valId-repmat(patt.^2, 1, 1, size(y, 3)/nbPatt)));
%% SIM Reconstruction
for ii=1:length(lamb_TV)
    FF=[Fn, {lamb_TV(ii)*F_TV}, {lamb_hess(ii)*F_hess}, {pos}];
    HH=[Hn, {Op_TV}, {Op_hess}, {OpD}];
    rho=[ones(size(Fn))*rhoDTNN,rho_TV,rho_hess,rhoDTNN];
    % Solver for ADMM
    A = rho_TV * Op_TV.makeHtH + rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(sz);
    solver = @(zn, rho_n, x0) solver_ADMM(A, HH, zn, rho_n);
    Opt=OptiADMM(F0, FF, HH, rho, solver);
    Opt.OutOp=MyOutputOpti(1,[],round(maxIt/10),1:(length(FF)-1));
    Opt.OutOp.saveXopt=0;
    Opt.CvOp=TestCvgStepRelative(1e-5);
    Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
    Opt.maxiter=maxIt;                % max number of iterations
    Opt.run(zeros(Hn{1}.sizein));     % run the algorithm
    %Opt.run(gt);     % run the algorithm
    if isGPU==1
        xopt=gather(Opt.xopt);
        fftxopt=gather(log(1+abs(fftshift(fftshift(Sfft(xopt,3),1),2))));   % because Sfft uses zeros_
    else
        xopt=Opt.xopt;
        fftxopt=log(1+abs(fftshift(fftshift(Sfft(xopt,3),1),2)));
    end
    if sav
        if ~isempty(pattname)
            saveastiff(single(xopt),[outFolder,...
                '_lamb_TV',num2str(lamb_TV(ii)),'_lamb_hess',num2str(lamb_hess(ii)),'.tif']);
            saveastiff(single(fftxopt),[outFolder,...
                '_lamb_TV',num2str(lamb_TV(ii)),'_lamb_hess',num2str(lamb_hess(ii)),'-FFT.tif']);
        else
            saveastiff(single(xopt),[outFolder,...
                '_lamb_TV',num2str(lamb_TV(ii)),'_lamb_hess',num2str(lamb_hess(ii)),'.tif']);
            saveastiff(single(fftxopt),[outFolder,...
                '_lamb_TV',num2str(lamb_TV(ii)),'_lamb_hess',num2str(lamb_hess(ii)),'-FFT.tif']);
        end
    end
end

% Rescale reconstruction
y_lambda = zeros(size(y));
for t = 1 : size(y, 3)
    y_lambda(:,:,t) = (S*H) * (squeeze(xopt(:,:,t)) .* patt(:,:,1+mod(t-1, size(patt, 3))));  
end
resc = sum(y_lambda(:) .* y(:)) / sum(y_lambda(:).^2);
y_rec = y_lambda * resc; im_rec = xopt * resc;

disp('###################################');
figure;sliceViewer(im_rec, 'Colormap', parula); title('Reconstructed image'); colorbar;
figure;sliceViewer(log(1+abs(fftshift(fftn(im_rec)))), 'Colormap', parula); title('Reconstructed image FFT');

%% Comparisons

if do_noTV
    % No TV regularization (frame by frame)
    Op_hess=LinOpHess(szUp,'circular');  % Hessian-Shatten: Hessian Operator
    F_hess=CostMixNormSchatt1(Op_hess.sizeout,1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
    pos=CostNonNeg(szUp);
    Fn={};Hn={};F0=[];
    xopt_noTV = zeros(sz);
    for t = 1 : size(y,3)
        L2=CostL2([],y(:,:,t));             % L2 cost function
        Fn=(L2*(S*H));
        Fn.doPrecomputation=true;
        Hn=LinOpDiag(H.sizein,patt(:,:,1+mod(t-1, size(patt, 3))));
        OpD=LinOpDiag(szUp,sqrt(valId-patt(:,:,1+mod(t-1, size(patt, 3))).^2));
        for ii=1:length(lamb_hess)
            FF=[{Fn}, {lamb_hess(ii)*F_hess}, {pos}];
            HH=[{Hn}, {Op_hess}, {OpD}];
            rho=[rhoDTNN,rho_hess,rhoDTNN];
            % Solver for ADMM
            A = rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(szUp);
            solver = @(zn, rho_n, x0) solver_ADMM(A, HH, zn, rho_n);
            Opt=OptiADMM(F0, FF, HH, rho, solver);
            Opt.OutOp=MyOutputOpti(1,[],round(maxIt/10),1:(length(FF)-1));
            Opt.OutOp.saveXopt=0;
            Opt.CvOp=TestCvgStepRelative(1e-5);
            Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
            Opt.maxiter=maxIt;                % max number of iterations
            Opt.run(zeros(Hn.sizein));     % run the algorithm
            xopt_noTV(:,:,t) = Opt.xopt;
        end
    end

    figure;sliceViewer(xopt_noTV, 'Colormap', parula); title('Reconstructed image (without TV)'); colorbar; caxis([min(y(:)), max(y(:))]);
    figure;sliceViewer(log(1+abs(fftshift(fftn(xopt_noTV)))), 'Colormap', parula); title('Reconstructed image (without TV) FFT');

end
