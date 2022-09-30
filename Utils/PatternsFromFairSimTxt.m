function [patt,orr,phOff,k0]=PatternsFromFairSimTxt(file, varargin)
% Generate pattern from FairSim estimated parameters.

%% Read file
data = fscanf(fopen(file), '%f %f %f', [3, Inf]);

%% Extract parameters
res_xy=64;   % Lateral resolution (converted to nm)
sz=512;          % Image size (square)
nbOrr=3;               % Number of orientations
nbPh=3;                % Number of phases
nbPatt=nbOrr*nbPh;  % Total number of patterns
nbTime = size(data, 2) / 3; % Number of time frames
phShift=linspace(0,2*pi,nbPh+1); phShift=phShift(1:end-1);                     % Phase shift
phOff = zeros(3, nbTime);
shift = zeros(2, 3, nbTime);
orr = zeros(3, nbTime);
for t = 1 : nbTime
    for ii=1:nbOrr
        phOff(ii, t) = data(3, ii + (t-1)*nbPh);
        shift(:, ii, t) = data(1:2, ii + (t-1)*nbPh);
        nrm = norm(shift(:, ii, t));
        sinA=shift(2, ii, t)/nrm;
        cosA=shift(1, ii, t)/nrm;
        orr(ii, t)=-sign(sinA)*acos(cosA);
    end
end

%% Take mean over time frame of all parameters

phOff = mean(phOff, 2); % Big variability
shift = mean(shift, 3);
orr = mean(orr, 2);
for ii=1:nbOrr
    k0(ii) = norm(shift(:, ii));
end

if nargin<=1
    szFinal=[sz,sz];
else
    szFinal=varargin{1};
end

%% Build patterns

patt=zeros([szFinal(1)*2,szFinal(2)*2,nbPatt]);
[X,Y]=meshgrid(0:szFinal(2)*2-1,0:szFinal(1)*2-1);X=X*res_xy/2;Y=Y*res_xy/2;
it=1;
for ii=1:length(orr)
    k=2*pi*[cos(orr(ii)), sin(orr(ii))]*(k0(ii))/(sz*res_xy*2);
    for jj=1:length(phShift)
        patt(:,:,it)=1+ cos(2*(k(1)*X+k(2)*Y + phOff(ii)/2+ phShift(jj)/2));
        it=it+1;
    end
end
% patt=patt-min(patt(:))+1e-3;
nbPatt=size(patt,3); % Normalization such that the mean of each pattern is 1/#Patterns
for ii=1:nbPatt
    tmp=patt(:,:,ii);
    patt(:,:,ii)=patt(:,:,ii)/(mean(tmp(:))*nbPatt);
end


end
