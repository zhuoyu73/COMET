function [ shot_navs ] = calcMRENavPhaseErrors_alt(rInfo, navigatorImages)
   %calcMRENavPhaseErrors References phase errors of navigators in preparation
   %for a pcSENSE reconstruction.
   %
   %  [ shot_navs ] = calcMRENavPhaseErrors(rInfo, navigatorImages);
   %
   %   Inputs:
   %     rInfo - RecoInfo object that has been inialized from the raw data for
   %             the reconstruction. We require this to have access to the number of
   %             images, phases, etc.
   %     navigatorImages  - Complex navigator images
   %
   %   Outputs:
   %     shot_navs - Referenced phase error maps for each shot at the resolution
   %                 of the imaging data.
   %
   
   shot_navs = zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes,rInfo.nRepetitions);
   img_nav_HR = zeros(size(shot_navs));
   
   for ii = 1:rInfo.nShots
       for jj = 1:rInfo.nPartitions
           for kk = 1:rInfo.nSlices
               for ll = 1:rInfo.nAverages
                   for mm = 1:rInfo.nPhases
                       for nn = 1:rInfo.nEchoes
                           for oo = 1:rInfo.nRepetitions
                               % modifying this to match how Joe had it
                               nav_data = ifftshift(ifftshift(fft2(fftshift(fftshift(navigatorImages(:,:,:,ii,jj,kk,ll,mm,nn,oo),1),2)),1),2);
                               nav_data_HR = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor);
                               nav_data_HR(floor((rInfo.N-rInfo.NNav)/2)+1:floor((rInfo.N-rInfo.NNav)/2)+rInfo.NNav,floor((rInfo.N-rInfo.NNav)/2)+1:floor((rInfo.N-rInfo.NNav)/2)+rInfo.NNav,:) = nav_data;
                               img_nav_HR(:,:,:,ii,jj,kk,ll,mm,nn,oo) = fftshift(fftshift(ifft2(ifftshift(ifftshift(nav_data_HR,1),2)),1),2);
                           end
                       end
                   end
               end
           end
       end
   end
   
   for ii = 1:rInfo.nShots
       for jj = 1:rInfo.nPartitions
           for kk = 1:rInfo.nSlices
               for ll = 1:rInfo.nAverages
                   for mm = 1:rInfo.nPhases
                       for nn = 1:rInfo.nEchoes
                           for oo = 1:rInfo.nRepetitions
                               %Calculate reference image
                               refImage = reshape(img_nav_HR(:,:,:,:,:,kk,ll,mm,nn,oo),rInfo.N,rInfo.N,rInfo.multibandFactor,[]);
                               refImage = mean(refImage,4);
                               shot_navs(:,:,:,ii,jj,kk,ll,mm,nn,oo) = angle(img_nav_HR(:,:,:,ii,jj,kk,ll,mm,nn,oo)./refImage);
                           end
                       end
                   end
               end
           end
       end
   end
   shot_navs(isnan(shot_navs))=0;
   
%    shot_navs = zeros(rInfo.nRepetitions,rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
%    parfor oo = 1:rInfo.nRepetitions
%       shot_navs_temp =  zeros(rInfo.N,rInfo.N, rInfo.multibandFactor, rInfo.nShots, rInfo.nPartitions, rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes);
%       for nn = 1:rInfo.nEchoes
%          for mm = 1:rInfo.nPhases
%             for ll = 1:rInfo.nAverages
%                for jj = 1:rInfo.nPartitions
%                   for kk = 1:rInfo.nSlices
%                      for ii = 1:rInfo.nShots
%                         for hh = 1:rInfo.multibandFactor
%                            shot_navs_temp(:,:,hh,ii,jj,kk,ll,mm,nn) = resample_map_resolution(shot_navsTemp(:,:,hh,ii,jj,kk,ll,mm,nn,oo),rInfo.N,1,rInfo.FOV*10,rInfo.FOV*10);
%                         end
%                      end
%                   end
%                end
%             end
%          end
%       end
%       
%       % Reassemble output variable
%       shot_navs(oo,:,:,:,:,:,:,:,:,:,:) = shot_navs_temp;
%    end
%    shot_navs = permute(shot_navs,[2:10,1]);
   
end
