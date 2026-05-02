function [imgPhantom, FM, FMdx, FMdy] = genTestCase(N, peakOffRes)
   %GENTESTCASE Generate the phantom for the test case
   %   Detailed explanation goes here
   
   DefineBrain;
   imgPhantom = RasterizePhantom(Brain,[N,N]);
   [FM, FMdx, FMdy] = genSynthFieldMap(N,peakOffRes);
   FM = FM.';
   FMdx = FMdx.';
   FMdy = FMdy.';
   
end

