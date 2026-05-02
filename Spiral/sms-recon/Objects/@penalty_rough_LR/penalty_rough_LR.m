function ob = penalty_rough_LR(Cn, v, ntp)
%
%
% Joe Holtrop
% University of Illinois
% October 4, 2013
  


%	default object
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);

if nargin == 0
	ob = class(ob, 'penalty_rough_LR');
	return
end

ob.v = v;
ob.ntp = ntp;
ob.R = size(v,1); %Rank of model being used
ob.C = Cn; %Rank of model being used

ob.is.empty	= logical(0);
ob = class(ob, 'penalty_rough_LR');






