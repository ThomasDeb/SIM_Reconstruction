clear; close all;
setGlobalBioImPath;

%% Simulation
% -- Paths
gtpath='./Data/object.tif';      % file name ground truth
gtsize = [256, 256];            % size of ground truth image
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

% -- Acquisition
photBud=0; %500;    % Photon Budget
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
sav=1;                                  % Boolean if true save the result

% -- Data
valback=0;       % Background value

% -- Objective Function
lamb_TV=1e-5;       % Hyperparameter (can be an array to loop)
lamb_hess=lamb_TV/3;       % Hyperparameter (can be an array to loop)
symPsf=1;        % Boolean true if psf is symmetric

% -- SIM reconstruction
alg=1;                    % Algorithm (1: ADMM / 2: Primal-Dual)
maxIt = 100;               % Max iterations
ItUpOut =round(maxIt/10);  % Iterations between to call to OutputOpti
rhoDTNN = 1e-3;            % rho parameter (ADMM) associated to data term
rho_TV = lamb_TV*10;             % rho parameter (ADMM) associated to TV term (must be greater or equal than rhoDTNN if iterCG=0)
rho_hess = lamb_hess*100;             % rho parameter (ADMM) associated to Hessian-Schatten term (must be greater or equal than rhoDTNN if iterCG=0)
split=2;                  % Splitting strategy for data fidelity:
valId=2;                  % Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 3)

%% Reconstruction
disp('########## Reconstruction #########');

global isGPU
addpath Utils

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
if alg==2,split=1;  end

%% Conversion CPU/GPU is necessary
psf=gpuCpuConverter(psf);
y=gpuCpuConverter(y);

%% Common operators and costs
% -- Data term
%psf=psf/sum(psf(:));
otf=fft2(fftshift(psf));
%H=LinOpConv(repmat(otf, [1, 1, nbPatt]), 1, 1:2);
H=LinOpConv(otf, 1, 1:2);
S=LinOpDownsample(szUp,downFact);
Fn={};Hn={};F0=[];
if ~isempty(pattname)
    nn=sum(patt.^2,3);                      % Trick to perform the last ADMM step in Fourier
    patt=patt/(sqrt(max(nn(:))));
    for i=1:size(patt,3)
        Si = LinOpSelectorPlanes(sz,3,i,true);
        L2=CostL2([],y(:,:,i));             % L2 cost function
        switch split
            case 0
                if i==1
                    F0=L2*S*H*LinOpDiag(H.sizein,patt(:,:,i))*Si;
                else
                    F0=F0+L2*S*H*LinOpDiag(H.sizein,patt(:,:,i))*Si;
                end
            case 1
                Fn{i}=L2;
                Fn{i}.doPrecomputation=true;
                Hn{i}=S*H*LinOpDiag(H.sizein,patt(:,:,i))*Si;
            case 2
                Fn{i}=(L2*(S*H));
                Fn{i}.doPrecomputation=true;
                Hn{i}=LinOpDiag(H.sizein,patt(:,:,i))*Si;
        end
    end
else
    L2=CostL2([],y);             % L2 cost function
    LS=(L2*(S*H));
    LS.doPrecomputation=true;
end
% -- Regularization
    % TV regularization in time
    Op_TV=LinOpGrad(sz,3,'circular');      % TV regularizer: Operator Gradient
    F_TV=CostL1(Op_TV.sizeout);                % TV regularizer: L1 norm (gradient of size 1)
    Op_hess=LinOpHess(sz,'circular',[1 2]);  % Hessian-Shatten: Hessian Operator
    F_hess=CostMixNormSchatt1(Op_hess.sizeout,1);  % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
% -- Non-Negativity constraint
pos=CostNonNeg(sz);
if alg==1 && (~isempty(pattname) && split==2)
    assert(valId >= 1,'valId must be greater than 1');
    OpD=LinOpDiag(sz,sqrt(valId-patt.^2));
else
    OpD=LinOpIdentity(sz);
end
%% SIM Reconstruction 



for ii=1:length(lamb_TV)
    if alg==1     % ADMM
        FF=[Fn, {lamb_TV(ii)*F_TV}, {lamb_hess(ii)*F_hess}, {pos}];
        HH=[Hn, {Op_TV}, {Op_hess}, {OpD}];
        rho=[ones(size(Fn))*rhoDTNN,rho_TV,rho_hess,rhoDTNN];
        % Solver for ADMM
        A = rho_TV * Op_TV.makeHtH + rho_hess * Op_hess.makeHtH + rhoDTNN * valId * LinOpIdentity(sz);
        solver = @(zn, rho_n, x0) solver_ADMM(A, HH, zn, rho_n);
        Opt=OptiADMM(F0, FF, HH, rho, solver);
        Opt.OutOp=MyOutputOpti(1,gt,round(maxIt/10),1:(length(FF)-1));
    elseif alg==2 % Primal-Dual
        FF=[Fn,{lamb(ii)/(size(otf,3))*Freg}];
        HH=[Hn,{Opreg}];
        Opt=OptiPrimalDualCondat(F0,pos,FF,HH);
        Opt.tau=tau;
        Opt.rho=rho;
        T=HH{1}.makeHtH(); for ll=2:length(HH), T=T+HH{ll}.makeHtH(); end
        Opt.sig=1/(tau*T.norm);
        Opt.OutOp=MyOutputOpti(1,gt,round(maxIt/10),2:(2+nbPatt));
    end
    Opt.OutOp.saveXopt=0;
    Opt.CvOp=TestCvgStepRelative(1e-5);
    if alg==1 && (split~=2 && iterCG~=0)
        Opt.CG.maxiter=iterCG;
    end
    Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
    Opt.maxiter=maxIt;                % max number of iterations
    Opt.run(zeros(Hn{1}.sizein));     % run the algorithm 
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

disp('###################################');
figure;sliceViewer(gt, 'Colormap', parula); title('Ground truth');
figure;sliceViewer(xopt, 'Colormap', parula); title('Reconstructed image');
figure;sliceViewer(log(1+abs(fftshift(fftn(xopt)))), 'Colormap', parula); title('Reconstructed image FFT'); 