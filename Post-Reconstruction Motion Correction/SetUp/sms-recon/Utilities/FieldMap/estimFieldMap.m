function [FM,FMsmoothed] = estimFieldMap(FMImages,mask,TEs,varargin)
%ESTIMFIELDMAP Summary of this function goes here
%   Detailed explanation goes here

%% Deal with optional input arguments
%p = inputParser();

%p.parse(varargin{:});


%% Calculate unwrapped phase images
% Bring last dimension to first
FMImages = shiftdim(FMImages,ndims(FMImages)-1);
sizePhaseImages = size(FMImages);
phaseImageDims = ndims(FMImages);
unwrappedPhaseImages = zeros(size(FMImages));
preludeOptions = struct('mask',mask);
mkdir unwrap_niis
preludeOptions.path2UnwrappedPhase = 'unwrap_niis/unwrappedPhase.nii';
for ii = 1:sizePhaseImages(1)
    %unwrappedPhaseImages(ii,:) = col(lapunwrap(reshape(phaseImages(ii,:).*col(mask).',sizePhaseImages(2:end))));
    %unwrappedPhaseImages(ii,:) = col(reshape(phaseImages(ii,:).*col(mask).',sizePhaseImages(2:end)));
    imageTemp = reshape(FMImages(ii,:),sizePhaseImages(2:end));
    unwrappedPhaseImages(ii,:) = col(prelude(angle(imageTemp),abs(imageTemp),preludeOptions));
end

%% Calculate intial fit
%FM = (backgroundField(end,:) - backgroundField(2,:))./(TEs(end)-TEs(2));
%dims = 1:ndims(backgroundField);
%backgroundField = permute(backgroundField, [dims(2:end),dims(1)]);
%backgroundField = reshape(backgroundField, [sizePhaseImages(2:end),sizePhaseImages(1)]);
FMinit = (unwrappedPhaseImages(end,:) - unwrappedPhaseImages(1,:))./(TEs(end)-TEs(1));
FMinit = reshape(FMinit,sizePhaseImages(2:end));

unwrappedPhaseImages = permute(unwrappedPhaseImages,[2:sizePhaseImages(end),1]);

%% Now we can fit field maps with minimal hassle, right?
R = Robj(logical(ones(size(mask))),'edge_type','tight','order',2,'beta',2^-3,'type_denom','matlab','potential','quad');
for ii = 1:sizePhaseImages(1)
    Gtemp = vertcat(spdiag(TEs(ii).*ones(prod(sizePhaseImages(2:end)),1),'nowarn'));
    if ii == 1
        G = Gtemp;
    else
        G = vertcat(G,Gtemp);
    end
end

FM = solve_pwls_pcg(col(FMinit),G,1,col(unwrappedPhaseImages),R,'niter',10);
FM = reshape(FM,sizePhaseImages(2:end));
if nargout > 1
    %% Expand using smoothn
    FM(~mask) = NaN;
    FMsmoothed = smoothn(FM);
    FM(isnan(FM)) = 0;
end
