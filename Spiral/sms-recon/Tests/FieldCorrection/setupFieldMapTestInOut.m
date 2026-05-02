%% Generate phantom and trajectory
clear
N = 128;
maxOffResonance = -200; %Hz
halfVoxelShiftOn = false;
timeSegSize = 1E-3; % Maximum Time Segment length in (s).
%Generate oversampled phantom and field map
OSFactor = 2;
[imgPhantom, FM, FMdx, FMdy] = genTestCase(OSFactor*N,maxOffResonance);

% Generate 18 shot spiral trajectory for field map calculation
nShots = 18;
maxGrad = 22; % mT/m
maxSlew = 120; %mT/m/s
FOV = 240;

[kx, ky]=genkspace(FOV*.1, N, 0, nShots, maxGrad*.1, maxSlew, 5E-6);
TEs = [0, 1E-3];
timingVec = (0:(length(kx)/nShots-1))*5E-6;
imgRange = [0,4.5];
%% Generate Reference Data

% Simulate asymmetric spin echo acquistion using timing vector

for ii = 1:length(TEs)
   fieldMapTimingVec(:,ii) = timingVec + TEs(ii);
   G = Gdft(ky,kx,zeros(size(kx)),OSFactor*N,OSFactor*N,1,2*pi*FM,squeeze(fieldMapTimingVec(:,ii)));
   rawData(:,ii) = G*col(imgPhantom);
   if halfVoxelShiftOn
      % Half Voxel Shift the data for reconstruction
      halfVoxelShift = -1/(2*N);
      rawData(:,ii) = rawData(:,ii).*exp(1j*(halfVoxelShift.*kx(:)+halfVoxelShift.*ky(:))*2*pi);
   end
end


save rawData.mat rawData imgPhantom FM N maxOffResonance nShots maxGrad maxSlew FOV
FMRecon = resample_map_resolution(FM,N,1,FOV,FOV);
FMReconDx = resample_map_resolution(FMdx,N,1,FOV,FOV);
FMReconDy = resample_map_resolution(FMdy,N,1,FOV,FOV);
%FMRecon = FM(1:OSFactor:end,1:OSFactor:end);
R = Robj(true(N,N),'edge_type','tight','order',2,'beta',0,'type_denom','matlab','potential','quad');
xinit = zeros(N,N);

[x, y] = meshgrid(1:N,1:N);
maskCirc = sqrt((x-(N+1)/2).^2+(y-(N+1)/2).^2)<((N)/2);


%% No Field Corrected NUFFT Recon and Field Map Estimation
reconCase = 'NUFFT Only Recon';

G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);

%A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,zeros(size(FMRecon)),timingVec,-1*Gy,-1*Gx,zeros(size(Gx)));
%A = Gdft(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,timingVec);
for ii = 1:length(TEs)
   imgEst(:,ii) = solve_pwls_pcg(col(xinit),G, 1, col(rawData(:,ii)), R, 'niter',10);
end

imgNUFFT = reshape(imgEst,N,N,length(TEs));
figure; im(imgNUFFT,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

% Generate mask based on the last echo
maskSmoothing = (col(abs(imgNUFFT(:,:,end))) > (0.1*max(col(abs(imgNUFFT(:,:,end))))));
maskSmoothing = reshape(maskSmoothing,N,N,1);
mask = maskSmoothing.*maskCirc;

fmNUFFT = -mri_field_map_reg(imgNUFFT,TEs,'l2b',-3,'mask',(mask>0));

%% Field Corrected Gdft Recon and Field Map Estimation
reconCase = 'Gdft Field Map Images';
G = Gdft(ky,kx,zeros(size(kx)),N,N,1,2*pi*FMRecon,squeeze(fieldMapTimingVec(:,1)));

for ii = 1:length(TEs)
   imgEst(:,ii) = solve_pwls_pcg(col(xinit),G, 1, col(rawData(:,ii)), R, 'niter',10);
end

imgGdft = reshape(imgEst,N,N,length(TEs));
figure; im(imgGdft,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

% Generate mask based on the last echo
maskSmoothing = (col(abs(imgGdft(:,:,end))) > (0.1*max(col(abs(imgGdft(:,:,end))))));
maskSmoothing = reshape(maskSmoothing,N,N,1);
mask = maskSmoothing.*maskCirc;

fmGdft = -mri_field_map_reg(imgGdft, TEs, 'l2b', -3, 'mask', (mask>0),'niter',1000);

%% Field Corrected Gdft_R2 Recon and Field Map Estimation
reconCase = 'Gdft with Gradients Field Map Images';
FMReconTest = fraccircshift(FMRecon,[1/(2*OSFactor),1/(2*OSFactor)]);
%Centered Differences in center, forward difference on edges
[Gx,Gy] = gradient(FMReconTest);
G = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,2*pi*FMRecon,zeros(size(FMRecon)),squeeze(fieldMapTimingVec(:,1)),Gx,Gy,zeros(size(Gx)));

for ii = 1:length(TEs)
   imgEst(:,ii) = solve_pwls_pcg(col(xinit),G, 1, col(rawData(:,ii)), R, 'niter',10);
end

imgGdftGrads = reshape(imgEst,N,N,length(TEs));
figure; im(imgGdftGrads,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

% Generate mask based on the last echo
maskSmoothing = (col(abs(imgGdft(:,:,end))) > (0.1*max(col(abs(imgGdft(:,:,end))))));
maskSmoothing = reshape(maskSmoothing,N,N,1);
mask = maskSmoothing.*maskCirc;

fmGdftGrads = -mri_field_map_reg(imgGdft, TEs, 'l2b', -3, 'mask', (mask>0),'niter',1000);
% %% Testing joint estimation of field map and images
% reconCase = 'Joint Estimation of Field Map and Images';
% 
% imgGdft = reshape(imgEst,N,N,length(TEs));
% figure; im(imgGdft,imgRange); title(reconCase);
% saveas(gcf,[reconCase '.png'], 'png');
% 
% % Generate mask based on the last echo
% maskSmoothing = (col(abs(imgGdft(:,:,end))) > (0.1*max(col(abs(imgGdft(:,:,end))))));
% maskSmoothing = reshape(maskSmoothing,N,N,1);
% mask = maskSmoothing.*maskCirc;
% L = ceil(max(timingVec)/timeSegSize);
% tt = repmat(fieldMapTimingVec,[nShots,1]);
% FMestim = FM_jointestimation(rawData,fmNUFFT,N,kx,ky,L,tt,[0,1E-3],2,15,10);
%% Field Corrected NUFFT Recon and Field Map Estimation

L = ceil(max(timingVec)/timeSegSize);
reconCase = ['Time Segmentation with L = ' num2str(L)];
G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
A = TimeSegmentation(G,squeeze(fieldMapTimingVec(:,1)),2*pi*FMRecon,L);

for ii = 1:length(TEs)
   imgEst(:,ii) = solve_pwls_pcg(col(xinit),A, 1, col(rawData(:,ii)), R, 'niter',10);
end

imgTimeSeg = reshape(imgEst,N,N,length(TEs));
figure; im(imgTimeSeg,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

% Generate mask based on the last echo
maskSmoothing = (col(abs(imgNUFFT(:,:,end))) > (0.1*max(col(abs(imgNUFFT(:,:,end))))));
maskSmoothing = reshape(maskSmoothing,N,N,1);
mask = maskSmoothing.*maskCirc;

fmTimeSeg = -mri_field_map_reg(imgTimeSeg,TEs,'l2b',-3,'mask',(mask>0),'niter',1000);

%% Let's estimate single shot images

% Generate 1 shot spiral trajectory for field map testing
nShots = 1;
maxGrad = 22; % mT/m
maxSlew = 120; %mT/m/s
FOV = 240;

[kx, ky]=genkspace(FOV*.1, N, 0, nShots, maxGrad*.1, maxSlew, 5E-6);
timingVec = (0:(length(kx)/nShots-1))*5E-6;

G = Gdft(ky,kx,zeros(size(kx)),OSFactor*N,OSFactor*N,1,2*pi*FM,timingVec(:),'VoxelBasis','boxcar');
rawData = G*col(imgPhantom);

if halfVoxelShiftOn
   % Half Voxel Shift the data for reconstruction
   halfVoxelShift = -1/(2*N);
   rawData = rawData.*exp(1j*(halfVoxelShift.*kx(:)+halfVoxelShift.*ky(:))*2*pi);
end

%% No Field Corrected NUFFT Recon
reconCase = 'NUFFT Only Recon';

G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
%A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,zeros(size(FMRecon)),timingVec,-1*Gy,-1*Gx,zeros(size(Gx)));
%A = Gdft(ky,kx,zeros(size(kx)),N,N,1,-2*pi*FMRecon,timingVec);

imgEst = solve_pwls_pcg(col(xinit),G, 1, rawData(:), R, 'niter',10);

imgNUFFTSingleShot = reshape(imgEst,N,N);
figure; im(imgNUFFTSingleShot,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

%% Field Corrected NUFFT Recon
L = ceil(max(timingVec)/timeSegSize);
reconCase = ['Time Segmentation with L = ' num2str(L)];
FMReconTest = fraccircshift(fmGdft,[-1/(2*OSFactor),-1/(2*OSFactor)]);

G = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
A = TimeSegmentation(G,timingVec,FMReconTest,L);
imgEst = solve_pwls_pcg(col(xinit),A, 1, rawData(:), R, 'niter',10);

imgTimeSegSingleShot = reshape(imgEst,N,N);
figure; im(imgTimeSegSingleShot,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

%% Gdft without Gradients Model Recon
reconCase = 'Gdft without Gradients Model';
FMReconTest = fraccircshift(fmGdft,[-1/(2*OSFactor),-1/(2*OSFactor)]);

A = Gdft(ky,kx,zeros(size(kx)),N,N,1,FMReconTest,timingVec);
imgEst = solve_pwls_pcg(col(xinit), A, 1, rawData(:), R, 'niter',10);

imgGdftSingleShot = reshape(imgEst,N,N);
figure; im(imgGdftSingleShot,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');

%% Gdft with Gradients Model Recon
reconCase = 'Gdft with Gradients Model';

%FMReconTest = resample_map_resolution(circshift(2*pi*FM,[1,1]),N,1,FOV,FOV);
%FMReconTest = fmGdft;
FMReconTest = fraccircshift(fmGdft,[-1/(2*OSFactor),-1/(2*OSFactor)]);
%Centered Differences in center, forward difference on edges
[Gx,Gy] = gradient(FMReconTest/(2*pi));
% Forward difference, don't care about edges for this test
%Gx = circshift(FMReconTest,[1,0]) - FMReconTest;
%Gy = circshift(FMReconTest,[0,1]) - FMReconTest;
%Gx = FMReconDy;
%Gy = FMReconDx;

A = Gdft_r2(ky,kx,zeros(size(kx)),N,N,1,FMReconTest,zeros(size(FMRecon)),timingVec,Gx,Gy,zeros(size(Gx)));
imgEst = solve_pwls_pcg(col(xinit), A, 1, rawData(:), R, 'niter',10);

imgGdftGradsSingleShot = reshape(imgEst,N,N);
figure; im(imgGdftGradsSingleShot,imgRange); title(reconCase);
saveas(gcf,[reconCase '.png'], 'png');