function [curve, roiImage] = checkExpFitROI(image,bValues)
figure;
im(image(:,:,1));

roi = imellipse;
sizeImages = size(image);
mask = roi.createMask();
meanValues = zeros(length(bValues),1);
for ii = 1:length(bValues)
    temp = image(:,:,ii);
    
    meanValues(ii) = mean(temp(mask));
end

[curve, goodness] = fit(bValues(:),meanValues,'exp1');
% Display fit stats
curve
% Display goodness
goodness

frame = getframe(gca);
roiImage = permute(rgb2gray(frame2im(frame)),[2,1]);

figure; plot(curve,bValues(:),meanValues);