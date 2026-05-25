function [r2,we] = createR2StarFM(rInfo, FM, sen, mask, varargin)
%createR2StarFM - Creates an R2* and Field Map from multiple echo
%                   time data.
%
% Syntax:  [R2Star, FM] = createR2StarFM(rInfo, FM, sen, mask, varargin);
%
% Inputs:
%    rInfo - recoInfo object that has been initialized with a Siemens
%    VB/VD/VE .dat file that can be read via mapVBVD
%    FM - An initial estimate of the field map calculated with
%    createFieldMap()
%    sen - SENSE Map of the same size as the FM estimate and images
%    mask - mask corresponding to one from SENSE map generation
%
% Options (Name/Value Pairs):
%    TEMask - A vector of logicals specifying whether or not the echo time
%    should be used. Defaults to all echoes.
%    UseHanning - Use Hanning window (true) or Min-Max window for time
%    segementation in field correction. Defaults to Hanning window (true).
% Outputs:
%    R2Star - Jointly estimated R2Star map of size [N,N,NSlices]
%    FM - Jointly estimated Field Map [N,N,NSlices]
%
% Example:
%    rInfo          = recoInfo(filename);
%    images         = gridCoilImages(rInfo);
%    [sen, mask]    = createSenMap(images,1);
%    [FM, FMImages] = createFieldMap(rInfo,sen,mask,1);
%    [R2Star, FM]       = createR2StarFM(rInfo, FM, sen, mask);
%
% Other m-files required: recoInfo.m, solve_pwls_pcg.m, Robject.m,
%       jointestimation_r2_we_multi_echo.m
% Subfunctions: none
% MAT-files required:
%
% Author: Giang-Chau Ngo
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% November 2016; Last revision: 5-Nov-2016

%% Deal with input parsing using inputParser Class
p = inputParser;

TEMask = logical(ones(rInfo.nEchoes,1));
UseHanning = true;
addParameter(p,'TEMask',TEMask);
addParameter(p,'UseHanning',UseHanning);
addParameter(p,'R2Init',[]);
addParameter(p,'InnerIterations',2);
addParameter(p,'OuterIterations',70);
p.parse(varargin{:});

TEMask = p.Results.TEMask;
UseHanning = p.Results.UseHanning;
R2Init = p.Results.R2Init;
InnerIterations = p.Results.InnerIterations;
OuterIterations = p.Results.OuterIterations;

%% Prep data and k-space trajectory

%Convert TEMask to indexes, hack to fix the fact that dataRead is a
%function pretending to be an array.
TEIndicies = 1:length(TEMask);
%Only handle 2D acquired case for now
data = rInfo.dataRead([],1,[],1,1,TEIndicies(TEMask),1,1);

%datatp = permute(data,[1 3 2 4]);
%use_hanning = 1;

for nsl = 1:rInfo.nSlices
%for nsl = 30
    if isempty(R2Init)
        % Not masked, just uses mask to fill in full size of image.
        r2init = double(mask(:,:,rInfo.ExcOrder(nsl))).*25; % inital value is 25Hz
    else
        r2init = R2Init(:,:,rInfo.ExcOrder(nsl));
    end
    [r2{nsl}, we{nsl}] = jointestimation_r2_we_multi_echo_fmincon(data(:,:,1,nsl,1,1,:,1,1,:)*3e3, TEMask, rInfo, FM(:,:,rInfo.ExcOrder(nsl)), sen(:,:,rInfo.ExcOrder(nsl),:), r2init, UseHanning, InnerIterations,OuterIterations,mask(:,:,rInfo.ExcOrder(nsl)));
    
end


end
