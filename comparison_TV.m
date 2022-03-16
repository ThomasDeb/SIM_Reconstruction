clear; close all;
setGlobalBioImPath;
addpath('./Utils/');
rng(1);

%% Simulation
% -- Paths
gtpath='./Data/object.tif';      % file name ground truth
gtsize = [512, 512];            % size of ground truth image
expFolder='./SimulatedExample/';   % experiment folder
sav=1;                           % to save results

% -- PSF
lamb=488;                % Illumination wavelength
res=32;                  % Resolution (nm)
Na=1.4;                  % Objective numerica aperture
nl=1.518;                % Refractive index of the objective medium (glass/oil)
ns=1.333;                % Refractive index of the sample medium (water)

% -- Patterns
orr=[0 pi/3 2*pi/3];   % Patterns orientations (vector)
ph=linspace(0,2*pi,4); % Patterns lateral phases (vector)
ph=ph(1:end-1);
a=0.9;                 % Amplitude coefficient
bet=asin(Na/nl);       % Angle between side beams and the optic axis (e.g. bet asin(Na/nl))
num_cycles = 3;

% -- Acquisition
photBud= 500; % 500;    % Photon Budget
downFact = [2 2];

% -- Run Simulator
disp('######## Generate SIM data ########');
run './SimuUFS.m'
disp('###################################');

%% Reconstruction
% -- Paths
basedir='./SimulatedExample/';
dataname=[basedir,'AcqData.tif'];      % File name data image
psfname=[basedir,'PSF'];               % File name psf
pattname=[basedir,'patterns.tif'];     % File name patterns
gtname=['./SimulatedExample/objectUFS.tif'];           % File name ground truth (if no let empty)
outFolder=[basedir];                    % Folder to save results
sav = true;                                  % Boolean if true save the result
do_noTV = false;
% -- Data
valback=0;       % Background value

% -- Objective Function
lamb_TV=5*1e-2;       % Hyperparameter (can be an array to loop)
lamb_hess=1*1e-2;       % Hyperparameter (can be an array to loop)
symPsf=1;        % Boolean true if psf is symmetric

% -- SIM reconstruction
maxIt = 100;               % Max iterations
ItUpOut = round(maxIt/10);  % Iterations between to call to OutputOpti
rhoDTNN = 1e-3;            % rho parameter (ADMM) associated to data term
rho_TV = lamb_TV/1000;             % rho parameter (ADMM) associated to TV term (must be greater or equal than rhoDTNN if iterCG=0)
rho_hess = lamb_hess/100;             % rho parameter (ADMM) associated to Hessian-Schatten term (must be greater or equal than rhoDTNN if iterCG=0)
valId = 2;                  % Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 3)

%% Reconstruction
disp('########## Reconstruction #########');

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
y=double(loadtiff(dataname))-valback; %maxy=max(y(:));y=y/maxy;              % Load and normalize data
% -- Ground Truth
if ~isempty(gtname)                                                        % Load ground truth if available
    gt=double(loadtiff(gtname))-valback;
else
    gt=[];
end
% -- Some parameters
szUp=size(psf);
sz=[szUp, nbPatt]; %input size

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
    for i=1:size(patt,3)
        Si = LinOpSelectorPlanes(sz,3,i,true);
        L2=CostL2([],y(:,:,i));             % L2 cost function
        Fn{i}=(L2*(S*H));
        Fn{i}.doPrecomputation=true;
        Hn{i}=LinOpDiag(H.sizein,patt(:,:,i))*Si;
    end
else
    L2=CostL2([],y);             % L2 cost function
    LS=L2*(S*H);
    LS.doPrecomputation=true;
end
% -- Regularization
% TV regularization in time
Op_TV_per = LinOpGrad(sz, 3, 'circular');      % TV regularizer: Operator Gradient
F_TV_per = CostL1(Op_TV_per.sizeout);                % TV regularizer: L1 norm (gradient of size 1)
Op_TV = LinOpIdentity(sz);
F_TV = CostTV1D(sz, 3);
Op_hess=LinOpHess(sz,'circular',[1 2]);  % Hessian-Shatten: Hessian Operator
F_hess=CostMixNormSchatt1(Op_hess.sizeout,1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
% -- Non-Negativity constraint
pos=CostNonNeg(sz);
assert(valId >= 1,'valId must be greater than 1');
OpD=LinOpDiag(sz,sqrt(valId-patt.^2));
%% SIM Reconstruction

for ii=1:length(lamb_TV)
    FF=[Fn, {lamb_TV(ii)*F_TV}, {lamb_hess(ii)*F_hess}, {pos}];
    HH=[Hn, {Op_TV}, {Op_hess}, {OpD}];
    rho=[ones(size(Fn))*rhoDTNN,rho_TV,rho_hess,rhoDTNN];
    % Solver for ADMM
    A = rho_TV * Op_TV.makeHtH + rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(sz);
    solver = @(zn, rho_n, x0) solver_ADMM(A, HH, zn, rho_n);
    Opt=OptiADMM(F0, FF, HH, rho, solver);
    Opt.OutOp=MyOutputOpti(1,im,round(maxIt/10),1:(length(FF)-1));
    Opt.OutOp.saveXopt=0;
    Opt.CvOp=TestCvgStepRelative(1e-5);
    Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
    Opt.maxiter=maxIt;                % max number of iterations
    Opt.run(zeros(Hn{1}.sizein));     % run the algorithm
    xopt=Opt.xopt;
    fftxopt=log(1+abs(fftshift(fftshift(Sfft(xopt,3),1),2)));
end

% Periodic TV
for ii=1:length(lamb_TV)
    FF_per=[Fn, {lamb_TV(ii)*F_TV_per}, {lamb_hess(ii)*F_hess}, {pos}];
    HH_per=[Hn, {Op_TV_per}, {Op_hess}, {OpD}];
    rho_per=[ones(size(Fn))*rhoDTNN,rho_TV,rho_hess,rhoDTNN];
    % Solver for ADMM
    A_per = rho_TV * Op_TV_per.makeHtH + rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(sz);
    solver_per = @(zn, rho_n, x0) solver_ADMM(A_per, HH_per, zn, rho_n);
    Opt_per=OptiADMM(F0, FF_per, HH_per, rho_per, solver_per);
    Opt_per.OutOp=MyOutputOpti(1,im,round(maxIt/10),1:(length(FF)-1));
    Opt_per.OutOp.saveXopt=0;
    Opt_per.CvOp=TestCvgStepRelative(1e-5);
    Opt_per.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
    Opt_per.maxiter=maxIt;                % max number of iterations
    Opt_per.run(zeros(Hn{1}.sizein));     % run the algorithm
    xopt_per=Opt_per.xopt;
    fftxopt_per=log(1+abs(fftshift(fftshift(Sfft(xopt_per,3),1),2)));
end



disp('###################################');
figure; sliceViewer(im, 'Colormap', parula); title('Ground truth'); colorbar;
figure; sliceViewer(xopt, 'Colormap', parula); title('Reconstructed image'); colorbar; caxis([min(im(:)), max(im(:))]);
figure; sliceViewer(xopt_per, 'Colormap', parula); title('Reconstructed image (periodic TV)'); colorbar; caxis([min(im(:)), max(im(:))]);
figure; sliceViewer(xopt - xopt_per, 'Colormap', parula); title('Difference (non)periodic TV'); colorbar;

SNR_UFS = OptSNR(im, xopt);
disp(['SNR = ', num2str(SNR_UFS),' dB']);

SNR_UFS_per = OptSNR(im, xopt_per);
disp(['SNR (periodic TV) = ', num2str(SNR_UFS_per),' dB']);

%% Comparisons

% Widefield image
SNR_widefield = OptSNR(im, repmat(mean(im, 3),[1, 1, nbPatt]));
disp(['SNR (noiseless widefield image) = ', num2str(SNR_widefield),' dB']);

if do_noTV
    % No TV regularization (frame by frame)
    Op_hess=LinOpHess(szUp,'circular');  % Hessian-Shatten: Hessian Operator
    F_hess=CostMixNormSchatt1(Op_hess.sizeout,1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
    pos=CostNonNeg(szUp);
    Fn={};Hn={};F0=[];
    xopt_noTV = zeros(sz);
    for i = 1 : size(patt,3)
        L2=CostL2([],y(:,:,i));             % L2 cost function
        Fn=(L2*(S*H));
        Fn.doPrecomputation=true;
        Hn=LinOpDiag(H.sizein,patt(:,:,i));
        OpD=LinOpDiag(szUp,sqrt(valId-patt(:,:,i).^2));
        for ii=1:length(lamb_hess)
            FF=[{Fn}, {lamb_hess(ii)*F_hess}, {pos}];
            HH=[{Hn}, {Op_hess}, {OpD}];
            rho=[rhoDTNN,rho_hess,rhoDTNN];
            % Solver for ADMM
            A = rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(szUp);
            solver = @(zn, rho_n, x0) solver_ADMM(A, HH, zn, rho_n);
            Opt=OptiADMM(F0, FF, HH, rho, solver);
            Opt.OutOp=MyOutputOpti(1,im(:,:,i),round(maxIt/10),1:(length(FF)-1));
            Opt.OutOp.saveXopt=0;
            Opt.CvOp=TestCvgStepRelative(1e-5);
            Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
            Opt.maxiter=maxIt;                % max number of iterations
            Opt.run(zeros(Hn.sizein));     % run the algorithm
            xopt_noTV(:,:,i) = Opt.xopt;
        end
    end

    figure;sliceViewer(xopt_noTV, 'Colormap', parula); title('Reconstructed image (without TV)'); colorbar; caxis([min(im(:)), max(im(:))]);
    figure;sliceViewer(log(1+abs(fftshift(fftn(xopt_noTV)))), 'Colormap', parula); title('Reconstructed image (without TV) FFT');

    SNR_noTV = OptSNR(im, xopt_noTV);
    disp(['SNR (without TV) = ', num2str(SNR_noTV),' dB']);
end
