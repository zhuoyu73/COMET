function[ mask, f, skull ] = bet(mag, varargin)
%BET - brain extraction tool
%
%   calls on BET of FSL package 
%       
%   Syntax
%
%   BET(Mag)
%   BET(Mag, Optional)
%
%   Description
%       (see BET man page)
%
%   Mask            = BET(Mag)                   
%       returns binary brain mask based on (single-echo) 3D magnitude input 
%
%   Mask, f         = BET(Mag, Options)
%       also returns additional fractional intensity threshold 
%
%   Mask, f, skull  = BET(Mag, Options)
%       also returns binary mask of skull (might not be useful)
%
%   The following Name-value pairs are supported
%
%       voxelSize
%               default: [1 1 1] (isotropic)
%
%       f 
%           "fractional intensity threshold"
%               default: 0.5
% 
%   2013 topfer@ualberta.ca 
% Found from: https://github.com/sunhongfu/QSM/
% Modified by Alex Cerjanic 10/16/2017 to change options handling to use
% name valuie pairs


%% constants

DEFAULT_VOXELSIZE = [2 2 2] ;
DEFAULT_F         = 0.3 ;

%% check inputs

if nargin < 1 || isempty(mag)
    error('Function requires at least 1 input.')
end
fValidationFcn = @(x) isscalar(x) && isnumeric(x) && (x <= 1) && (x>0);
voxelSizeValidationFcn = @(x) ~isscalar(x) && (ndims(x) <= 3) && (ndims(x)>1);
p = inputParser();

addOptional(p,'f', DEFAULT_F,fValidationFcn);
addOptional(p,'VoxelSize', DEFAULT_VOXELSIZE, voxelSizeValidationFcn);

parse(p,varargin{:});

voxelSize = p.Results.VoxelSize;
f = p.Results.f;

if nargin < 2
    disp('Default parameters will be used')
    disp('assuming isotropic resolution')
end

tmpFldr = 'tmpBETFldr/' ;
system(['mkdir ' tmpFldr]) ;

save_nii( make_nii( mag, voxelSize ), [tmpFldr 'mag.nii'] ) ;

A = mag;
tot_mass = sum(A(:));
[xx,yy,zz] = ndgrid(1:size(A,1),1:size(A,2),1:size(A,3));
RR = sum(xx(:).*A(:))/tot_mass;
CC = sum(yy(:).*A(:))/tot_mass;
SS = sum(zz(:).*A(:))/tot_mass;

betArgStr = [tmpFldr 'masked.nii -m -g 0.1 -w 1.1 -c ' num2str(round(RR)+15) ' '  num2str(round(CC)) ' ' num2str(round(SS)) ' -f ' num2str(f) ] ;

if nargout == 3
    betArgStr = strcat( betArgStr, ' -s ' ) ;
end
    
system(['$FSLDIR/bin/bet2 ' tmpFldr 'mag.nii ' betArgStr] ) ;

system(['mv ' tmpFldr 'masked.nii_mask.nii.gz ' tmpFldr 'masked_mask.nii.gz']) ; 
system(['gunzip ' tmpFldr 'masked_mask.nii.gz -d']) ;
mask = load_nii([tmpFldr 'masked_mask.nii']) ;
mask = single( mask.img ) ;

if nargout == 3
    system(['gunzip ' tmpFldr 'masked_skull.nii.gz -d']) ;
    skull = load_nii([tmpFldr 'masked_skull.nii']) ;
    skull = single( skull.img ) ;
end


system(['rm -r ' tmpFldr]) ;



end