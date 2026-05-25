function [SNR, gfactor, pseudoreplicas] = calcPseudoreplicaSNR(rInfo, sen, mask, varargin)
   % calcPseudoreplicaSNR Calculate SNR and gfactor via the pseudareplica
   % method. Note that field correction is not supported.
   %
   % Syntax:  [img] = calcPseudoreplicaSNR(rInfo,sen,mask, optional args)
   %
   % Inputs:
   %    rInfo   - intialized recoInfo object, rInfo must have the coil noise
   %              covariance matrix calcualted and stored in the recoInfo
   %              object.
   %    sen     - SENSE Maps corresponding to images to be reconstructed
   %              If the SENSE maps are not of the same size, they need to be
   %              interpolated before calling this function.
   %    mask    - Masks corresponding to images to be reconstructed
   %              If the masks are not of the same size, they need to be
   %              interpolated before calling this function.
   %
   % Optional Name - Value Pair Inputs:
   %    Rbeta   - Add beta for regularization function. Default Beta is zero
   %              which is not recommended.
   %
   %    Penalty - String recognized by Robject (or Robj) to specify penalty
   %              weighting. Default is 'quad'. Options, such as 'huber'
   %              exist.
   %
   %    delta   - Some penalties require a delta option, such as 'huber'. This
   %              has no effect if the penalty specified does not use the
   %              delta value.
   %
   %    Niter   - Sets the number of CG iterations used to solve for the image
   %              by solve_pwls_pcg.m. Default is 10.
   %
   %    NReplicas - Number of pseudoreplicas to compute. Will received N+1
   %                images back, as the first image in the pseudoreplicas
   %                output will be the reference image with no synthetic noise
   %                added. Default is 50.
   %   slicesToRecon - Array of indicies that correspond to images to be
   %                   reconstructed. Indicies can be in any order or
   %                   discontinuous. Default is 1:rInfo.nSlices.
   %
   %   averagesToRecon - Array of indicies that correspond to images to be
   %                     reconstructed. Indicies can be in any order or
   %                     discontinuous. Default is 1:rInfo.nAverages.
   %
   %   phasesToRecon - Array of indicies that correspond to images to be
   %                   reconstructed. Indicies can be in any order or
   %                   discontinuous. Default is 1:rInfo.nPhases.
   %
   %   echoesToRecon - Array of indicies that correspond to images to be
   %                   reconstructed. Indicies can be in any order or
   %                   discontinuous. Default is 1:rInfo.nEchoes.
   %
   %   repetitionsToRecon - Array of indicies corresponding to images to be
   %                        reconstructed. Indicies can be in any order or
   %                        discontinuous. Default is 1:rInfo.nRepetitions.
   %
   % Outputs:
   %    SNR    - SNR map calculated via pseudoreplica method
   %
   %    gfactor - gfactor map calculated via pseudoreplica method
   %
   %    pseudoreplicas - stack of NReplicas+1 images. (Pseudoreplicas(1) is the
   %                     reference image with no noise added.)
   %
   % Example:
   %    [cImages, sen, FM, FMImages] = ReconSenFM();
   %    rInfo  = rInfo(filename);
   %    images = gridCoilImages(rInfo);
   %    im(sum(abs(images).^2,10)); % Sum of Squares image for
   %                                            % first echo time in field map
   %
   %
   % Other m-files required: Too numerous to list
   % Subfunctions: none
   % MAT-files required: none
   %
   % Author: Alex Cerjanic
   % University of Illinois at Urbana-Champaign
   % email address:
   % Website:
   % 12-Sep-2017; Last revision: 12-Sep-2017
   
   %% Deal with Optional Inputs
   
   % Establish Defaults
   
   Rbeta = 0; %Use a zero beta, not a good idea in practice.
   RBetaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
   
   delta = 0; % Use a zero delta, not a good idea for huber.
   DeltaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
   
   penalty = 'quad'; % Specify defalt penalty weighting
   PenaltyValidationFcn = @(x) isstring(x);
   
   Niter = 10;
   NIterValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   
   dims2penalize = [1 1 0];
   dims2penalizeValidationFcn = @(x) isnumeric(x);
   
   slicesToRecon = 1:rInfo.nSlices;
   slicesToReconValidationFcn = @(x) isnumeric(x);
   
   averagesToRecon = 1:rInfo.nAverages;
   averagesToReconValidationFcn = @(x) isnumeric(x);
   
   phasesToRecon = 1:rInfo.nPhases;
   phasesToReconValidationFcn = @(x) isnumeric(x);
   
   echoesToRecon = 1:rInfo.nEchoes;
   echoesToReconValidationFcn = @(x) isnumeric(x);
   
   repetitionsToRecon = 1:rInfo.nRepetitions;
   repetitionsToReconValidationFcn = @(x) isnumeric(x);
   
   p = inputParser();
   addOptional(p,'Rbeta',Rbeta,RBetaValidationFcn);
   addOptional(p,'delta',delta,DeltaValidationFcn);
   addOptional(p,'penalty',penalty,PenaltyValidationFcn);
   addOptional(p,'L',L,LValidationFcn);
   addOptional(p,'Niter', Niter, NIterValidationFcn);
   addOptional(p,'dims2penalize', dims2penalize, dims2penalizeValidationFcn);
   addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
   addOptional(p,'averagesToRecon', averagesToRecon, averagesToReconValidationFcn);
   addOptional(p,'phasesToRecond', phasesToRecon, phasesToReconValidationFcn);
   addOptional(p,'echoesToRecon', echoesToRecon, echoesToReconValidationFcn);
   addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);
   
   parse(p,varargin{:});
   
   % Override default parameters with their successfully parsed values.
   Rbeta = p.Results.Rbeta;
   delta = p.Results.delta;
   penalty = p.Results.penalty;
   L = p.Results.L;
   Niter = p.Results.Niter;
   dims2penalize = p.Results.dims2penalize;
   
   %% Setup Recon
   
   % Note that we can perform the recon over an arbitrary set of images,
   % echoes, etc. To support this, we index over a vector of indicies. By
   % default the vector of inidices 1:nIndex, but could be something like
   % [2,4,5], etc if the user supplies those indexes.
   pseudoreplicas = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,...
      length(slicesToRecon),length(averagesToRecon),length(phasesToRecon),...
      length(echoesToRecon),length(repetitionsToRecon),nReplicas);
   
   for kk = 1:length(slicesToRecon)
      for ll = 1:length(averagesToRecon)
         for mm = 1:length(phasesToRecon)
            for nn = 1:length(echoesToRecon)
               for oo = 1:length(repetitionsToRecon)
                  
                  % Now we need to look up the image indexes
                  slc = slicesToRecon(kk);
                  avg   = averagesToRecon(ll);
                  phs   = phasesToRecon(mm);
                  eco   = echoesToRecon(nn);
                  rep   = repetitionsToRecon(oo);
                  
                  if rInfo.multibandFactor > 1
                     G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        col(rInfo.ky(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        col(rInfo.kz(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        rInfo.N,...
                        rInfo.N,...
                        rInfo.multibandFactor,...
                        'mask', logical(mask(:,:,:,slc)));
                     
                  else % 2D Case
                     G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        col(rInfo.ky(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        col(rInfo.kz(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        rInfo.N, ...
                        rInfo.N,...
                        1,...
                        'mask',logical(mask(:,:,:,slc)));
                  end
                  
                  R = Robj(logical(squeeze(mask(:,:,:,slc))), ...
                     'edge_type','tight', ...
                     'order',2,...
                     'beta',Rbeta,...
                     'type_denom','matlab',...
                     'potential',penalty,...
                     'delta',delta,...
                     'dims2penalize',dims2penalize);
                  
                  sen_tmp = reshape(sen(:,:,:,slc,:),rInfo.N*rInfo.N*rInfo.multibandFactor,[]);
                  S = sense(G,sen_tmp(logical(col(mask(:,:,:,slc))),:));
                  
                  % Calc Lower L matrix from choleski decomposition
                  % of the noise covariance matrix
                  if rInfo.preWhiten == 1
                     Lmatrix = 1; % Noise covariance is unity by definition
                  else
                     Lmatrix = chol(rInfo.noiseCorr,'lower');
                  end
                  
                  for pp = 1:NReplicas
                     %% Calculating a pseudoreplica
                     data = rInfo.dataRead([],[],kk,ll,mm,nn,oo);
                     
                     % Add noise for pseudoreplica
                     
                     if pp ~= 1
                        data = reshape(data,[],rInfo.nCoils).';
                        data = data + Lmatrix*randn(size(data));
                        data = data.';
                     end
                     tic
                     
                     % Setup Iterative Reconstruction
                     
                     data = col(data);
                     
                     % Setup initial estimate of all zeros.
                     imginit = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor);
                     xinit = col(imginit(logical(squeeze(mask(:,:,:,slc)))));
                     
                     % Solve for image
                     colReplica = solve_pwls_pcg(xinit, S, 1, data, R,...
                        'niter',Niter);
                     % deal with support mask to create square image
                     % Note that we are using the loop indicies since we
                     % need indicies that go from 1:Nimages, etc.
                     pseudoreplicas(:,:,:,kk,ll,mm,nn,oo) = ...
                        embed(colReplica, logical(mask(:,:,:,slc)));
                     
                     toc
                  end
               end
            end
         end
      end
   end
   
   % To Do: Calculate mean of pseudoreplicas - the first (no noise added) image
   % to get SNR images. Calculate the SD of the pseduoreplicas and divide by
   % sqrt(Rfactor) to get the gfactor.
   
end

