
%% Simulating Field Correction results
% This script will generate an analytical brain slice as well as an analytical
% field map in order to compare different methods of field correction. This
% method compares: no field correction, time segmentation, field corrected DSFT
% with zero order field map, and field corrected DSFT with first order field map
% model.

%% Generate phantom and trajectory
% Target resolution
N = 128;
halfVoxelShiftOn = false;
% Generate oversampled phantom and field map
OSFactor = 2;
timeSegSize = 1E-3; % Maximum Time Segment length in (s).
maxOffResonance = -200; %Hz
if 0
    load ImageTest.mat
    imgPhantom = ImageTest;
    load FMTest.mat
    FM = FMTest/(2*pi);
else
    [imgPhantom, FM] = genTestCase(OSFactor*N,maxOffResonance);
end
% Generate single shot spiral trajectory and encoding parameters 
nShots = 1;
maxGrad = 22; % mT/m
maxSlew = 120; %mT/m/s
FOV = 240;

[kx, ky]=genkspace(FOV*.1, N, 0, nShots, maxGrad*.1, maxSlew, 5E-6); 

timingVec = (0:(length(kx)/nShots-1))*5E-6;
if 0
    timingVec = zeros(size(timingVec));
end
imgRange = [0,4.5*max(abs(imgPhantom(:)))];
R = Robj(true(N,N),'edge_type','tight','order',2,'beta',1E-1,'type_denom','matlab','potential','quad');
xinit = zeros(N,N);

%% Generate Reference Data

G = Gdft(ky,kx,zeros(size(kx)),OSFactor*N,OSFactor*N,1,2*pi*FM,timingVec);

rawData = G*col(imgPhantom);

if halfVoxelShiftOn
   % Half Voxel Shift the data for reconstruction 
   halfVoxelShift = 1/(2*N);
   rawData = rawData.*exp(1j*(halfVoxelShift.*kx(:)+halfVoxelShift.*ky(:))*2*pi);
end
save rawData.mat rawData imgPhantom FM N maxOffResonance nShots maxGrad maxSlew FOV
FMRecon = resample_map_resolution(FM,N,1,FOV,FOV);
%FMRecon = wiener2(FMRecon,[3,3]);

%% No Field Corrected NUFFT Recon
reconCase = 'NUFFT Only Recon';

G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
%A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,zeros(size(FMRecon)),timingVec,-1*Gy,-1*Gx,zeros(size(Gx)));
%A = Gdft(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,timingVec);

[imgEst, ~, normsNUFFT] = solve_pwls_pcg(col(xinit),G, 1, rawData(:), R, 'niter',10,'isave','all');

imgNUFFT = reshape(imgEst,N,N,[]);
figure; im(imgNUFFT,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');


%% Field Corrected NUFFT Recon
L = ceil(max(timingVec)/timeSegSize);
reconCase = ['Time Segmentation with L = ' num2str(L)];
%FMReconTest = fraccircshift(FMRecon,[1/(2*OSFactor),1/(2*OSFactor)]);
FMReconTest = fraccircshift(FMRecon,[0,0]);
G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
A = TimeSegmentation(G,timingVec,2*pi*FMReconTest,L);
[imgEst,~,normsTimeSeg] = solve_pwls_pcg(col(xinit),A, 1, rawData(:), R, 'niter',10,'isave','all');

imgTimeSeg = reshape(imgEst,N,N,[]);
figure; im(imgTimeSeg,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

%% Gdft without Gradients Model Recon
reconCase = 'Gdft without Gradients Model';
%FMReconTest = fraccircshift(FMRecon,[1/(2*OSFactor),1/(2*OSFactor)]);
FMReconTest = fraccircshift(FMRecon,[0,0]);

A = Gdft(ky,kx,zeros(size(kx)),N,N,1,2*pi*FMReconTest,timingVec);
[imgEst,~,normsGdft] = solve_pwls_pcg(col(xinit),A, 1, rawData(:), R, 'niter',10,'isave','all');

imgGdft = reshape(imgEst,N,N,[]);
figure; im(imgGdft,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

%% Gdft with Gradients Model Recon
reconCase = 'Gdft with Gradients Model';
%FMReconTest = fraccircshift(FMRecon,[1/(2*OSFactor),1/(2*OSFactor)]);
FMReconTest = fraccircshift(FMRecon,[0,0]);
[Gx,Gy] = gradient(FMReconTest);
A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,2*pi*FMReconTest,zeros(size(FMRecon)),timingVec,Gy,Gx,zeros(size(Gx)));
%A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,2*pi*FM128,zeros(size(FMRecon)),timingVec,zeros(size(FM128Dy)),zeros(size(FM128Dx)),zeros(size(Gx)));
[imgEst,~,normsGdftGrads] = solve_pwls_pcg(col(xinit),A, 1, rawData(:), R, 'niter',10,'isave','all');

imgGdftGrads = reshape(imgEst,N,N,[]);
figure; im(imgGdftGrads,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

