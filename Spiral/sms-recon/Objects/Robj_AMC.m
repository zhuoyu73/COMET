 function R = Robj_AMC(kappa, varargin)
%function R = Robj(kappa, [options])
%
% Build roughness penalty regularization "object" based on C = Cdiff() object,
% for regularized solutions to inverse problems.
%
% General form of nonquadratic penalty function:
%	R(x) = \sumk w_k \pot([Cx]_k), where [Cx]_k = \sum_j c_{kj} x_j.
%	
% For quadratic case, \potk(t) = t^2/2 in which case
%	R(x) = x' C' W C x / 2, where W depends on beta and edge_type.
%
% Penalty gradient is C' D C x, where D = diag{\wpot_k([Cx]_k)}.
%
% in
%	kappa	[nx,ny[,nz]]	kappa array, or logical support mask
%
% options
%	edge_type '?'		how to handle mask edge conditions (see Cdiff)
%		'none'		no roughness penalty (NOT DONE)
%		'tight'		only penalize within-mask differences (default)
%		'leak'		penalize between mask and neighbors
%				(its primary use for consistency with ASPIRE)
%	'order', ?		1st-order or 2nd-order differences (see Cdiff)
%	'offsets', [?]		offsets to neighboring pixels
%					(see Cdiff for the defaults)
%				use '3d:26' to penalize all 13 pairs of nbrs
%				use '0' for C = I (identity matrix)
%	'beta', ?		global regularization parameter
%					default: 2^0
%	'delta', ?		potential parameter, see potential_func()
%					or {delta, param}.  default: inf
%	'potential', '?'	e.g., 'huber', see potential_func()
%				default: 'quad' for quadratic regularization.
%	'type_denom', '?'	type of "denominator"
%					(for quadratic surrogates like SPS)
%		'matlab'	denominator for SPS
%		'aspire'	denominator for SPS that matches aspire
%		'none'		no denominator (default)
%				(because R.E and R.denom needed only for SPS)
%	'distance_power', ?	1 classical (default), 2 possibly improved
%					see penalty_mex()
%	'user_wt', [?]		User-provided array of penalty weight values
%					size: [nx,ny[,nz],#offsets]
%				of dimension [size(mask) length(offsets)].
%				These are .* the usual wt values for edge_type.
%	'mask'			Override default: mask = (kappa ~= 0)
%
% out
%	R structure has the following inline "methods" 
%	R.penal(R, x)	evaluates R(x)
%	R.cgrad(R, x)	evaluates \cgrad R(x) (column gradient)
%	R.denom(R, x)	evaluates denominator for separable surrogate
%	[pderiv pcurv] = feval(R.dercurv, R, C1*x) derivatives and curvatures
%			for non-separable parabola surrogates
%	R.diag(R)	diagonal of Hessian of R (at x=0), for preconditioners.
%	R.C1		differencing matrix, with entries 1 and -1,
%			almost always should be used in conjunction with R.wt
%
% Typical use:	mask = true([128 128]); % or something more conformal
%		R = Robject(mask, 'beta', 2^7);
%
% Copyright 2004-11-14, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(kappa, 'test'), run_mfile_local 'Robject_test', return, end

% option defaults
R.potential = 'quad';
R.beta = 2^0;
R.delta = inf;
R.edge_type = 'tight';
R.type_denom = 'none';
R.distance_power = 1;
R.order = 1; % 1st-order differences, used in call to Cdiff
R.user_wt = []; % user-provided wt values
R.mask = [];
R.offsets = [];
R.dims2penalize = [];

% parse name/value option pairs
R = vararg_pair(R, varargin);

% dimensions, and default offsets
if isempty(R.offsets)
	if ndims(kappa) == 2
		[nx ny] = size(kappa);
		R.offsets = [1 nx nx+1 nx-1];
	elseif ndims(kappa) == 3
		[nx ny nz] = size(kappa);
		R.offsets = [1 nx nx*ny];
	else
		error 'only 2D and 3D done'
	end
elseif streq(R.offsets, '3d:26') % all 26 neighbors (13 pairs)
	[nx ny nz] = size(kappa);
	R.offsets = [1 nx+[0 1 -1] nx*ny+col(outer_sum([-1:1],[-1:1]*nx))'];
end

if ~iscell(R.delta)
	R.delta = {R.delta};
end

R.offsets = int32(R.offsets);

if streq(R.edge_type, 'none')
	error 'todo: unpenalized'
end

if R.beta < 0, warn('Negative beta? This is probably wrong!'), end

%R.isquad = streq(R.potential, 'quad');

if isempty(R.mask)
	R.mask = kappa ~= 0; % default is to infer from kappas
end
if ~islogical(R.mask), error 'mask must be logical', end

% incorporate user provided wk values
if ~isempty(R.dims2penalize)
	R.dims2penalize = 1:ndims(R.mask);
end

%
% build plain sparse differencing object "C1", containing 1 and -1 values
% (or identity matrix)
%

[R.C1, R.wt] = C_sparse(R.mask,R.dims2penalize);
R.wt = R.beta * R.wt(:); % absorb beta into wk factors

% incorporate user provided wk values
if ~isempty(R.user_wt)
	R.wt = R.wt .* R.user_wt(:);
	R.user_wt = []; % save memory
end

% trick: for quadratic case, provide a R.C object where C = sqrt(W) * C1
%if R.isquad
	%R.C = R.C1;
	%R.C.cascade_after = diag_sp(sqrt(R.wt));
%end


%
% desired potential function
%

R.pot.delta = R.delta{:};
if isvar('delta')
	pot.delta = delta(:);
	clear delta
end

% trick: huber2 is just a huber with delta / 2
% so that weighting function drops to 1/2 at delta, like hyper3 etc.
if streq(R.potential, 'huber2')
	type = 'huber';
	R.pot.delta = R.pot.delta / 2;
end

% trick: hyper3 is just a hyperbola with delta scaled by sqrt(3)
% to approximately "match" the properties of 'cauchy' (old erroneous 'hyper')
if streq(R.potential, 'hyper3')
	type = 'hyper2';
	R.pot.delta = R.pot.delta / sqrt(3);
end

switch R.potential

%
% quadratic potential function
%
case 'quad'
	potk = '(abs(t).^2)/2';
	wpot = 'ones(size(t))';
	dpot = 't';

%
% huber potential function
%
case 'huber'
	potk = 'huber_pot(t, pot.delta)';
	wpot = 'huber_wpot(t, pot.delta)';
	dpot = 'huber_dpot(t, pot.delta)';

%
% cauchy penalty: d^2 / 2 * log(1 + (t/d)^2)
% Not convex!
%
case 'cauchy'
	potk = 'pot.delta.^2 / 2 .* log(1 + abs(t ./ pot.delta).^2)';
	wpot = '1 ./ (1 + abs(t ./ pot.delta).^2)';
	dpot = 't ./ (1 + abs(t ./ pot.delta).^2)';

%
% Geman&McClure penalty: d^2 / 2 * (t/d)^2 / (1 + (t/d)^2)
% Not convex!
%
case 'geman&mcclure'
	potk = 'pot.delta.^2 / 2 .* (t/pot.delta)^2 ./ (1 + abs(t ./ pot.delta).^2)';
	wpot = '1 ./ (1 + abs(t ./ pot.delta).^2).^2';
	dpot = 't ./ (1 + abs(t ./ pot.delta).^2).^2';

%
% hyperbola penalty: d^2 * [ sqrt(1 + (t/d)^2) - 1 ]
%
case 'hyper2'
	potk = 'pot.delta.^2 .* (sqrt(1 + abs(t ./ pot.delta).^2) - 1)';
	wpot = '1 ./ sqrt(1 + abs(t ./ pot.delta).^2)';
	dpot = 't ./ sqrt(1 + abs(t ./ pot.delta).^2)';

case 'hyper'
	error 'use "cauchy" or "hyper3" not "hyper" now'

%
% Lange1 penalty
%
case 'lange1'
	potk = 't.^2 / 2 ./ (1+abs(t./pot.delta))';
	wpot = '(1 + abs(t ./ pot.delta) / 2) ./ (1 + abs(t ./ pot.delta)).^2';
	dpot = ['t .* (' wpot ')'];

%
% Lange3 penalty
%
case 'lange3'
	potk = 'pot.delta.^2 .* (abs(t./pot.delta) - log(1+abs(t./pot.delta)))';
	wpot = '1 ./ (1 + abs(t ./ pot.delta))';
	dpot = 't ./ (1 + abs(t ./ pot.delta))';

%
% li98cfs
%
case 'li98cfs'
	% f = inline('atan(x) / x - 0.5'); fsolve(f, 2.3)
	R.pot.delta = R.pot.delta / 2.3311;
	potk = 'pot.delta.^2 .* ((t ./ pot.delta) .* atan(t ./ pot.delta) - 0.5 * log(1 + (t ./ pot.delta).^2))';
	wpot = 'li98cfs_wpot(t, pot.delta)';
	dpot = ['t .* ' wpot];

otherwise
		fail('Unknown potential "%s"', type)
end

if ~isfield(R.pot, 'potk')
	R.pot.potk = inline(potk, 'pot', 't');
	R.pot.wpot = inline(wpot, 'pot', 't');
	R.pot.dpot = inline(dpot, 'pot', 't');
end

%
% inline functions
%
R.cgrad = @(R, x) R.C1' * (diag_sp(R.wt) * R.pot.dpot(R.pot, R.C1 * x));
R.penal = @(R, x) ... % 2014-04-27 added 'double' to better match new version
	sum(repmat(R.wt, ncol(x)) .* R.pot.potk(R.pot, R.C1 * x), 'double');
% R.cgrad = @(R, x) R.C1' * (sparse(double(R.wt),1:length(R.wt),1:length(R.wt)) * R.pot.dpot(R.pot, R.C1 * x));
R.denom = @(R, ddir, x) ((R.C1*ddir).*(R.wt>0))' * (R.pot.wpot(R.pot,(R.C1*x).*(R.wt>0)) .* ((R.C1*ddir).*(R.wt>0)));
