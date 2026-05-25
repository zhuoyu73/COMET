function [bvalues] = CalcBValuesFromXml(filename)
% Function to calculate b-values
% Alex Cerjanic - 06/04/2017


DiffusionInfo = xml2struct(filename);

numDiffusionEncodings = str2double(DiffusionInfo.diffusionScheme.diffEncodings.Text);
diffEncodingDirection = str2num(DiffusionInfo.diffusionScheme.diffEncodingDirection.Text);

for ii = 1:numDiffusionEncodings
   diffGradients(ii) = str2num(DiffusionInfo.diffusionScheme.diffEncodingTable.diffEncode{1,ii}.Attributes.G);
   diffBigDelta_us(ii) = str2num(DiffusionInfo.diffusionScheme.diffEncodingTable.diffEncode{1,ii}.Attributes.BigDelta);
   diffLittleDelta_us(ii) = str2num(DiffusionInfo.diffusionScheme.diffEncodingTable.diffEncode{1,ii}.Attributes.LittleDelta);
end

% Hack to Fix 

%if exist('
RiseTime = 1520E-6;
% Calculate b-values 

gamma = 2*pi*42.57747892*1E6; %Radians per Tesla for Proton

bVal = @(G,BigDelta,LittleDelta,RiseTime)( gamma.^2.*(G/1E3).^2.*(LittleDelta./1E3).^2 .* ((BigDelta./1E3 - LittleDelta./(3.*1E3)) + (RiseTime/1E3).^3) - (LittleDelta./1E3).*((RiseTime./1E3).^2)./6)./1E6; 

bvalues = bVal(diffGradients,diffBigDelta_us/1E3,diffLittleDelta_us/1E3,0.4);


