% initializeContribPaths.m - Alex Cerjanic 2017/10/15 
% Adds externally developed tools to MATLAB paths. When you add tools to the
% external folder, add the required paths to this file to ensure they are
% available to all users.
   
% Get the root path of the repository
rootPathContrib = fileparts(mfilename('fullpath'));

% Add NIfTI Toolbox
addpath([ rootPathContrib '/NIfTI']);

% Add EPFL Analytical Brain Phantom
addpath(genpath([ rootPathContrib '/Phantoms']));

% Add miscelaneous single files
addpath([ rootPathContrib '/Misc']);

% Add Laplace Boundary Value Code files
addpath([ rootPathContrib '/LBV']);

% Add Multidimensional DCT files
addpath([ rootPathContrib '/DCT']);

% Add 3D Density Compensation Function
addpath([ rootPathContrib '/sdc3_nrz_11aug']);

% Add Brian Hargreave's modified mintgrad code
addpath([ rootPathContrib '/mintgrad']);

% Add ADiGator for automatic differentiation
addpath([ rootPathContrib '/adigator']);
startupadigator

% Add natsort for nice sorting of strings
% Used for post PowerGrid file processing
addpath([ rootPathContrib '/natsortfiles']);

