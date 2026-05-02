function [SNRApproximate] = quantifySNR(testImage)

im(testImage);

disp('Choose region representing Signal');
roiSNR = imellipse(gca);
maskSNR = roiSNR.createMask();
disp('Choose region representing background');

roiBackground = imellipse(gca);
maskBackground = roiBackground.createMask();

meanSNR = mean(testImage(maskSNR));

stdBackground = std(testImage(maskBackground));

SNRApproximate = meanSNR/stdBackground*.655;

