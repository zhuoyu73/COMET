function [D,rsquared] = checkExpFit(image,bValues)

figure;
im(image(:,:,1));

roi = imellipse;
mask = roi.createMask().';

sizes = size(image(:,:,1));
D = zeros(sizes);
A = zeros(sizes);
rsquared = zeros(sizes);

for ii = 1:sizes(1)
    for jj = 1:sizes(2)
        if mask(ii,jj)
            [curve, goodness] = fit(bValues(:),squeeze(image(ii,jj,:)),'exp1');
            D(ii,jj) =  -1*curve.b;
            A(ii,jj) = curve.a;
            rsquared(ii,jj) = goodness.adjrsquare;
        end
    end
end
% Compute Stats about Fit in ROI
stdD = std(D(mask));
meanD = mean(D(mask));

stdD * 1E3
meanD * 1E3

meanA = mean(A(mask));

bPlot = linspace(min(bValues),max(bValues),100);
curveFit = meanA.*exp(-1*meanD.*bPlot);

bValsUnique = unique(bValues);
for ii = 1:length(unique(bValues))
    tempImage = image(:,:,bValues == bValsUnique(ii));
    imageBvalue = tempImage(repmat(mask,[1,1,size(tempImage,3)]));
    meanImage(ii) = mean(col(imageBvalue));
end
figure;
plot(bPlot,curveFit);
hold on
%dataToPlot = image(repmat(mask,1,1,length(bValues)));
%bValuesToPlot = repmat(permute(bValues,[2,3,1]),sum(mask(:)),1);
%plot(bValuesToPlot(:),dataToPlot,'b.');
plot(bValsUnique,meanImage,'b.');
[curve, goodness] = fit(bValsUnique(:),col(meanImage),'exp1');
curve
goodness

figure;
histogram(D(mask))