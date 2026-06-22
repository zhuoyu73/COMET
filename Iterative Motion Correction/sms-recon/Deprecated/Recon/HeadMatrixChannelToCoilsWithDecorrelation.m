function [dataCoils] = HeadMatrixChannelToCoilsWithDecorrelation(ascconv,data,NoiseScan,nl,nadc)

%% Figuring out the right order of the coils
reorderArray = zeros(1,12);

for ii = 1:12 %NCoils
    %Cluster Number (1:4 for 12Ch Head Matrix)
    ClusterNumber = str2num(ascconv.asCoilSelectMeas.asList(ii).sCoilElementID.tElement(2));
    switch ascconv.asCoilSelectMeas.asList(ii).sCoilElementID.tElement(3)
        case 'P'
            reorderArray((ClusterNumber-1)*3 + 1) = ii;
        case 'S'
            reorderArray((ClusterNumber-1)*3 + 2) = ii;
        case 'T'
            reorderArray((ClusterNumber-1)*3 + 3) = ii;
        otherwise
            error('CoilElementID in ascconv header does not match expected form. Are you using this function with a coil other than the 12 Ch Head Matrix?');
            keyboard
    end
    
end

%% Preparing decorrelating the noise using the noise scan
sizeNoiseScan = size(NoiseScan);
NoiseCorr = 1/(sizeNoiseScan(2)-1)*(NoiseScan*NoiseScan');

L = chol(NoiseCorr,'lower');
Linv = inv(L);

% Decorrelating the noise by multiplying by the inverse of the square root 
% of the noise correlations
dataSize = size(data);
data = permute(data,[4,1,3,2]);
for ii = 1:nl
    for jj = 1:nadc
        dataWhitened(:,:,jj,ii) = Linv*data(:,:,jj,ii);
    end
end

dataWhitened = permute(dataWhitened,[2,4,3,1]);
dataWhitened = dataWhitened(:,:,:,reorderArray);

%% Transforming from the Tim channels to coils
for ii = 1:4 %Number of clusters of 3 channels
    %Left Channel
    dataCoils(:,:,:,(ii-1)*3+1) = 0.5.*(dataWhitened(:,:,:,(ii-1)*3+1) + sqrt(2).*dataWhitened(:,:,:,(ii-1)*3+2) + dataWhitened(:,:,:,(ii-1)*3+3));
    %Middle Channel
    dataCoils(:,:,:,(ii-1)*3+2) = 1j.*1./sqrt(2).*(dataWhitened(:,:,:,(ii-1)*3+1) - dataWhitened(:,:,:,(ii-1)*3+3));
    %Right Channel
    dataCoils(:,:,:,(ii-1)*3+3) = 0.5.*(sqrt(2).*dataWhitened(:,:,:,(ii-1)*3+2) - dataWhitened(:,:,:,(ii-1)*3+1) - dataWhitened(:,:,:,(ii-1)*3+3));
end

%dataCoils = permute(dataCoils,[3,2,1,4]);
end