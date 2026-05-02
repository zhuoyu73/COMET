function ob = LRob(A, v, index_list)
%function ob = fast_mr(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo,epiflag)
%	Construct MRI object, which can do Ax and A'y operations
% typical input params:
%  A: a cell array of objects 
%  v: the temporal basis function [R ntp]
%  A_indexes [optional]: array of indexes to map each frame to an object
%
% Joe Holtrop
% University of Illinois
% October 1, 2013
  


%	default object
ob.A = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);

if nargin == 0
	ob = class(ob, 'LRobj');
	return
end

%	fill object

ob.A = A;
ob.v = v;
ob.ntp = size(v,2);
ob.R = size(v,1); %Rank of model being used
if nargin < 3
    index_list = 1:size(v,2);
end

%check vector sizes
if size(v,2) ~= length(index_list)
    keyboard
end

ob.A_index_list=index_list;

ob.is.empty	= logical(0);

ob = class(ob, 'LRobj');






