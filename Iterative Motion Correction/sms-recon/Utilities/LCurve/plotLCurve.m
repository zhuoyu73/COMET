function plotLCurve(resid,roughness,lambdas )
   %plotLCurve Summary of this function goes here
   %   Detailed explanation goes here
   
   
   plot(resid,roughness);
   hold on
   for ii = 1:length(lambdas)
      plot(resid(ii),roughness(ii),'ro');
      text(resid(ii),roughness(ii),['\lambda=' num2str(lambdas(ii))]);
   end
      
 
end

