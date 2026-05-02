 function ob = Gnufft(varargin)
%function ob = Gnufft([mask,] args)
%	Construct Gnufft object, which computes nonunform DSFT samples
%	of signals with dimensions [(Nd)] via the NUFFT.
%
%	The arguments are simply a cell array of the all the arguments
%	that will be passed to "nufft_init()" in the appropriate order.
%	See Gnufft_test.m for example usage.
%
%	Basically, you create a system matrix object by calling:
%		G = Gnufft( ... )
%	and then you can use it thereafter by typing commands like
%		y = G * x;
%	which will auto-magically evaluate the DSFT samples.
%	This is useful for iterative image reconstruction in MRI.
%
%	Besides simple utilities like display, there are the following
%	capabilities of this object:
%	y = G * x		forward operation
%	x = G' * y		adjoint operation
%	
%	Optional arguments for which a value must follow:
%		'mask'		binary support array
%
%	Copyright 2003-6-1	Jeff Fessler	The University of Michigan

%
%	default object
%
warning('Deprecated: This object is deprecated and not recommended for new MRI reconstruction routines. Use NUFFT instead.');

ob.dims = [0 0];
ob.apower	= 1;		% array power, becomes 2 for G.^2
%ob.scale	= 1;		% global scaling factor used in some children
ob.mask		= [];		% [(Nd)] logical array
ob.chat		= logical(0);
ob.index1	= [];		% row indices (unless transposed) e.g. G(3:9,:)
ob.index2	= [];		% col indices (unless transposed) e.g. G(:,3:9)
				% the empty default means *all* rows/cols
ob.is_transpose = logical(0);
ob.is_masked	= logical(0);	% set to 1 if G(:,mask(:))
ob.is_subref	= logical(0);	% set to 1 if subscripted

ob.st = [];
ob.Nd = [];

if nargin < 1
	warning 'Gnufft called with too few arguments'
	help Gnufft
	ob = class(ob, 'Gnufft');
	return

elseif nargin == 1
	ob.mask = [];	% logical(ones(nx,ny));
	args = varargin{1};

elseif nargin == 2
	ob.mask = varargin{1};
	if ~islogical(ob.mask), error 'mask must be logical', end
	args = varargin{2};

else
	error 'too many arguments'
end

if ~isa(args, 'cell'), error 'args must be cell array', end

omega = args{1};
if (length(args{2}) == 1)
  ob.Nd = [args{2} 1];
else
  ob.Nd = args{2};
end
ob.dims = [size(omega,1) prod(ob.Nd)];

%
%	initialize the structure needed by NUFFT
%
ob.st = nufft_init(args{:});

ob.is.empty = logical(0);
ob = class(ob, 'Gnufft');
