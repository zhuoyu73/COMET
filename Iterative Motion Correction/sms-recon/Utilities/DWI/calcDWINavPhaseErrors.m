function [ phase_errors, shot_navs, mag_errors ] = calcDWINavPhaseErrors(rInfo, navigatorImages, tform)
%calcDWINavPhaseErrors Prepares phase errors of navigators in preparation
%for a pcSENSE reconstruction.
%
%  [ shot_navs ] = calcDWINavPhaseErrors(rInfo, navigatorImages);
%
%   Inputs:
%     rInfo - RecoInfo object that has been inialized from the raw data for
%             the reconstruction. We require this to have access to the number of
%             images, phases, etc.
%     navigatorImages  - Complex navigator images
%
%   Outputs:
%     shot_navs - Phase error maps for each shot at the resolution
%                 of the imaging data.
%

%resizedShots = zeros(size(navigatorImages));
resizedShots = zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes,rInfo.nRepetitions);

%% Resample navigators to full resolution
shot_navs = zeros(rInfo.nRepetitions,rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
blurKern = 1;
for oo = 1:rInfo.nRepetitions
    for nn = 1:rInfo.nEchoes
        for mm = 1:rInfo.nPhases
            for ll = 1:rInfo.nAverages
                for jj = 1:rInfo.nPartitions
                    for kk = 1:rInfo.nSlices
                        for ii = 1:rInfo.nShots
                            %for hh = 1:rInfo.multibandFactor
                            temp = navigatorImages(:,:,:,ii,jj,kk,ll,mm,nn,oo);
                            %smoothedShots(:,:,hh,ii,jj,kk,ll,mm,nn) = imgaussfilt(real(temp),blurKern) + 1j*(imgaussfilt(imag(temp),blurKern));
                            %resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo) = resample_map_resolution(temp,rInfo.N,[],rInfo.FOV*10,rInfo.FOV*10);
                            resizedShots(:,:,:,ii,jj,kk,ll,mm,nn,oo) = resampleNavigators(temp,rInfo);
                            %end
                        end
                    end
                end
            end
        end
    end
end

%% Register imagesim(angle(testSmoothed))
% Register the B0 navigator to the b0 images if we present a tform cell
% array
if (nargin > 2)
    for oo = 1:rInfo.nRepetitions
        for nn = 1:rInfo.nEchoes
            for mm = 1:rInfo.nPhases
                for ll = 1:rInfo.nAverages
                    for kk = 1:rInfo.nSlices
                        for jj = 1:rInfo.nPartitions
                            for ii = 1:rInfo.nShots
                                for hh = 1:rInfo.multibandFactor
                                    tempTform = tform{hh,kk};
                                    tempReal = imwarp(real(resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo)),tempTform,'OutputView',imref2d([rInfo.N,rInfo.N]));
                                    tempImag = imwarp(imag(resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo)),tempTform,'OutputView',imref2d([rInfo.N,rInfo.N]));
                                    
                                    resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo) = tempReal + 1j*tempImag;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Perform Image smoothing
for oo = 1:rInfo.nRepetitions
    shot_navs_temp =  zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
    
    for nn = 1:rInfo.nEchoes
        for mm = 1:rInfo.nPhases
            for ll = 1:rInfo.nAverages
                for jj = 1:rInfo.nPartitions
                    for kk = 1:rInfo.nSlices
                        for ii = 1:rInfo.nShots
                            for hh = 1:rInfo.multibandFactor
                                %shot_navs_temp(:,:,hh,ii,jj,kk,ll,mm,nn) = resample_map_resolution(smoothedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo),rInfo.N,[],rInfo.FOV*10,rInfo.FOV*10);
                                temp = resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn,oo);
                                if(blurKern == 0)
                                    shot_navs_temp(:,:,hh,ii,jj,kk,ll,mm,nn,oo) = temp;
                                else
                                    shot_navs_temp(:,:,hh,ii,jj,kk,ll,mm,nn,oo) = imgaussfilt(real(temp),blurKern) + 1j*(imgaussfilt(imag(temp),blurKern));
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Reassemble output variable
    shot_navs(oo,:,:,:,:,:,:,:,:,:,:) = shot_navs_temp;
end
shot_navs = permute(shot_navs,[2:10,1]);

%% Reference Phase
phase_errors = zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes,rInfo.nRepetitions);

% Create B_0 reference images via averages
b0Images = reshape(shot_navs(:,:,:,:,:,:,1,1,1,1),rInfo.N,rInfo.N,rInfo.multibandFactor,[],rInfo.nSlices);
b0RefImages = squeeze(mean(b0Images,4)); % Create a complex average of all of the shots/phases

for oo = 1:rInfo.nRepetitions
    %shot_navs_temp =  zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
    
    for nn = 1:rInfo.nEchoes
        for mm = 1:rInfo.nPhases
            for ll = 1:rInfo.nAverages
                for jj = 1:rInfo.nPartitions
                    for kk = 1:rInfo.nSlices
                        for ii = 1:rInfo.nShots
                            %for hh = 1:rInfo.multibandFactor
                                %temp = resizedShots(:,:,hh,ii,jj,kk,ll,mm,nn);
                                %shot_navs_temp(:,:,hh,ii,jj,kk,ll,mm,nn) = prelude(angle(temp),maskMB(:,:,jj,ii));
                                %refImage = reshape(shot_navs(:,:,:,:,:,kk,ll,mm,nn,oo),rInfo.N,rInfo.N,rInfo.multibandFactor,[]);
                                %refImage = mean(refImage,4);
                                phase_errors(:,:,:,ii,jj,kk,ll,mm,nn,oo) = angle(shot_navs(:,:,:,ii,jj,kk,ll,mm,nn,oo)./b0RefImages(:,:,:,kk)); %AMC
                            %end
                        end
                    end
                end
            end
        end
    end
end

%% Reference Magnitude
mag_errors = zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes,rInfo.nRepetitions);

% Create B_0 reference images via averages
%b0Images = reshape(shot_navs(:,:,:,:,:,:,1,1,1,1),rInfo.N,rInfo.N,rInfo.multibandFactor,[],rInfo.nSlices);
%b0RefImages = squeeze(mean(abs(b0Images),4)); % Create a magnitude average of all of the shots/phases

for oo = 1:rInfo.nRepetitions
    %shot_navs_temp =  zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
    
    for nn = 1:rInfo.nEchoes
        for mm = 1:rInfo.nPhases
            for ll = 1:rInfo.nAverages
                for kk = 1:rInfo.nSlices
                    refImage = reshape(shot_navs(:,:,:,:,:,kk,ll,mm,nn,oo),rInfo.N,rInfo.N,rInfo.multibandFactor,[]);
                    refImage = mean(abs(refImage),4);
                    for jj = 1:rInfo.nPartitions
                        for ii = 1:rInfo.nShots
                            %for hh = 1:rInfo.multibandFactor
                            mag_errors(:,:,:,ii,jj,kk,ll,mm,nn,oo) = abs(shot_navs(:,:,:,ii,jj,kk,ll,mm,nn,oo))./refImage; %AMC
                            %end
                        end
                    end
                end
            end
        end
    end
end


%% prep output
% Reassemble output variable
%     shot_navs(oo,:,:,:,:,:,:,:,:,:,:) = shot_navs_temp;
% end
%shot_navs = permute(shot_navs,[2:10,1]);

%phase_errors = shot_navs;
%phase_errors = angle(shot_navs);
phase_errors(isnan(phase_errors))=0;
% Use a loss function to try and control magnitude weighting (both +/-).
mag_errors = 1-abs(tanh(mag_errors));
mag_errors(~isfinite(mag_errors))=1;
end
