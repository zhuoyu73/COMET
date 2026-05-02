function [ resid, roughness, images ] = calcLCurve(reconFcn,lambdaVec)
   %calcLCurve Calculates L curve for a given reconstruction problem.
   %   Detailed explanation goes here
   %   [ resid, roughness, images ] = calcLCurve(reconFcn,lambdaVec);
   % Inputs:
   %  reconFcn  - a Function Handle that takes only a beta value to run, see
   %              example, it must return [image, resid, roughness]
   %  lambdaVec - a vector of regularization (lambda) values to calculate
   %              the l-curve along
   %
   % Outputs:
   %  resid     - residuals (norm(y-A*x,2)) for the problem at each beta (lambda
   %  roughness - rougness penalties (norm(R*x,2)) for the final image for each
   %              lambda value
   %  images    - Images calculated for each lambda value for inspection
   %
   % Example:
   %
   %  reconFunction = @(lambda) fieldCorrectedReconMB(rInfo, senMB, maskMB, ...
   %                                   FMMB, 'Niter', 30, 'L', 0, 'Rbeta', ...
   %                                   lambda, 'slicesToRecon', 4, ...
   %                                   'repetitionsToRecon', 4);
   %  decadesLambda = -5:1:5;
   %  [ resid, roughness, images ] = calcLCurve(reconFunction, 10.^decadesLambda);
   %  plotLCurve(resid,roughness,10.^decadesLambda);
   
   resid = zeros(size(lambdaVec));
   roughness = zeros(size(lambdaVec));
   images = cell(size(lambdaVec));
   
   for ii = 1:length(lambdaVec)
      [images{ii}, roughnessTemp, resid(ii)] = reconFcn(lambdaVec(ii));
      roughness2 = roughnessTemp{end};
      roughness(ii) = roughness2(end);
   end
   
   % Note that Robject absorbs beta into the rougness operator. We need to
   % divide by the beta to ensure that we are only looking at roughness and not
   % beta*roughness.
   roughness = roughness./lambdaVec;
   
end

