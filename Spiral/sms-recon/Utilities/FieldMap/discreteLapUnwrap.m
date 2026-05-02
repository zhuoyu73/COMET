function [unwrappedPhaseImage] = discreteLapUnwrap(phaseImage)
%discreteLapUnwrap Implements FUDGE unwrapping of phase images
%   See: Ching Lih Ng A, Pyng Ng M, Josan S, Farquharson S, Mulcahy C, Ordidge RJ. "Fast Unwrapping using Discrete Gradient Evaluation (FUDGE): an
%       analytical correction to the Laplacian-based phase unwrapping
%       technique for discrete data." ISMRM Annual Symposium 2016, #0027
%
%   Input:
%       phaseImages - 2D or 3D phase image on [-pi,pi] to be unwrapped


% Step 1: Compute the discrete second order difference (laplacian) operator
% in each of the 2 or 3 directions
sizeImages = size(phaseImage);
laplaceDifferences = zeros(sizeImages);

for ii = 1:ndims(phaseImage)
    rotArray = zeros(1,ndims(phaseImage));
    rotArray(ii) = 1;
    forDiff = circshift(phaseImage,rotArray) - phaseImage;
    backDiff = phaseImage - circshift(phaseImage,-1*rotArray);
    
    laplaceDifferences = laplaceDifferences + atan(sin(forDiff)./cos(forDiff)) - atan(sin(backDiff)./cos(backDiff));
end    
    
    
% Step 2: Compute the unwrapped difference using the DCT (including the
% multidimensional DCT code from MATLAB File Exchange, found in our Contrib
% folder).

dctLapPhase = mirt_dctn(laplaceDifferences);
[X,Y,Z] = meshgrid((1:sizeImages(1)),(1:sizeImages(2)),(1:sizeImages(3)));
cosTerm = (2*cos(pi/sizeImages(1).*(X))-2) + (2*cos(pi/sizeImages(2).*Y)-2) + (2*cos(pi/sizeImages(3).*Z)-2);
unwrappedPhaseImage = mirt_idctn(1./(cosTerm).*dctLapPhase);

end

