classdef pcSENSE_3D_Zhuoyu
    %pcSENSE Phase Corrected SENSE object for use correcting motion induced
    %phase errors (Version 2.0)
    %
    % Syntax:  S = pcSENSE(A,phaseMaps, senMap, varargin)
    %
    % Inputs:
    %    A       - TimeSegmentation (or compatible) object. For the nested
    %              G object, k-space should be stored as a 2D array with
    %              dimensionality [NreadoutPts, Nshots].
    %       
    %    phaseMaps - Real valued phase maps with range: [-2npi:pi] of size
    %                [Nx,Ny,Nz,nShot]
    %    senMaps   - Complex valued sensitivity maps with size
    %                [Nx,Ny,Nz,nCoils]
    %   
    % Outputs:
    %    S     - pcSENSE class object with overloaded mtimes operator for
    %            computing the forward and adjoint transform operations through
    %            A*x and A'*x.
    %
    % Example:
    %    % Create a spin-warp kspace
    %    [kx,ky] = ndgrid(-63:64,-63:64);
    %    % Initialize the object for a 2D 128x128 matrix recon.
    %    G =
    %    NUFFT(kx(:),ky(:),zeros(size(kx(:))),128,128,1);
    %    % Compute a data domain estimate
    %    imageModel = phantom(128,128);
    %    dataEstimate = G*imageModel(:);
    %    imageEstimate = G'*dataEstimate; 
    %    im(reshape(imageEstimate,128,128));
    %
    %
    % Other m-files required: nufft implementation and private directories
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Alex Cerjanic (Version 2.0)
    %         Ahn Van (Version 1.0)
    %         Joseph Holtrop (Version 1.0)
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website: www.mrfil.bioen.illinois.edu
    % 8-Sep-2017; Last revision: 8-Sep-2017    
    properties
        isEmpty = true;
        is2DPhaseMaps = false; % Are the phase maps 2D (for thin, multislab case)?
        sen = []; % Sense map
        phaseMaps = []; % Phase (real valued [-n*pi,n*pi]) maps for correction of motion induced phase error
        % Phase maps must be same size as imaging data with
        % exception of the 3d dimension (partition
        % dimension). If partition dimension is 1, then
        % Phase Maps will be repmat'ed up to the size of
        % the imaging volume in that dimension only.
        isTranspose = false; % Required to implement forward/adjoint transforms via mtimes
        A = 1; % Transform object supporting mtimes overloads
        nShots; % Number of shots to be corrected
        A_Array; % Array of shotwise transform objects to be created.
        Nimage; %Size of the image as a column vector
        ShotLength; %Size of k-space for a single shot as a column vector.
        nCoils; % Number of coils
        %Zhuoyu add input T
        T = [];
    end
    
    methods
        %Zhuoyu add input T
        function obj = pcSENSE_3D_Zhuoyu(A,phaseMaps,sen,T,varargin)
        %function obj = pcSENSE(A,phaseMaps,sen,varargin)
            
            % Check to make sure that the transform object has kspace
            % in the correct format.
   %         if (size(A.G.Kx,2) == 1)
   %             error('The A object passed in must have kspace stored in a 2D array [npts x nshots].');
   %         end
            
            obj.nShots = size(A.G.Kx,2);
            obj.ShotLength = size(A.G.Kx,1);
            obj.sen = sen;
            obj.nCoils = size(obj.sen,2);
            obj.phaseMaps = phaseMaps;
            %Zhuoyu add input T
            obj.T = T;

            if size(phaseMaps,3) == 1
                obj.is2DPhaseMaps = true;
            end
            obj.A = A;
            if obj.A.G.isMasked
                obj.Nimage = sum(obj.A.G.mask(:));
            else
                obj.Nimage = (obj.A.G.Nx*obj.A.G.Ny*obj.A.G.Nz);
            end
            
            obj.A.G.st = []; % Intentionally clear out the spatial interpolator
            % Intialize the field corrected transform objects for shot-wise
            % reconstruction
            for ii = 1:obj.nShots
                % Deal with support mask on the transform
                if A.G.isMasked
                    Gtemp = NUFFT(obj.A.G.Kx(:,ii),obj.A.G.Ky(:,ii),obj.A.G.Kz(:,ii),obj.A.G.Nx,obj.A.G.Ny,obj.A.G.Nz,'InterpMethod',obj.A.G.InterpMethod ,'mask',obj.A.G.mask);
                else % No support mask used
                    Gtemp = NUFFT(obj.A.G.Kx(:,ii),obj.A.G.Ky(:,ii),obj.A.G.Kz(:,ii),obj.A.G.Nx,obj.A.G.Ny,obj.A.G.Nz,'InterpMethod',obj.A.G.InterpMethod);
                end
                % Deal with Timing vector
                if (length(obj.A.timingVec) ~= length(obj.A.G.Kx(:,1)))
                    timingVec = obj.A.timingVec(length(obj.A.G.Kx(:,1)));
                else
                    timingVec = obj.A.timingVec;
                end
                obj.A_Array{ii} = TimeSegmentation(Gtemp,timingVec,obj.A.fieldMap,0);
                
                % Avoid the cost of recomputing the temporal interpolators,
                % let's just use them as is from the prototype object.
                obj.A_Array{ii}.timeInterp = obj.A.timeInterp;
                obj.A_Array{ii}.nShots = 1; % By definition since we are correcting phase, shot-wise
                obj.A_Array{ii}.L = obj.A.L;
                obj.A_Array{ii}.Interpolator = obj.A.Interpolator;
                
            end
            
            obj.isEmpty = false;
            
        end
        
        function y = size(obj)
            y = size(obj.A.G);
            y(1) = obj.nCoils*y(1);
        end
        
        function obj = ctranspose(obj)
            obj.isTranspose = ~obj.isTranspose;
        end
        
        function y = mtimes(obj,x)
            
            if ~obj.isTranspose % Forward operation
                y = zeros(obj.ShotLength*obj.nCoils, obj.nShots);
                %x_3d = reshape(x, [120, 120, 4]); % ZS reshape x to 120x120x4
                x_embed = embed(x, obj.A.G.mask);
                disp(size(x_embed)); %120x120x4

                for indx = 1:obj.nShots
                    Psense = zeros(obj.Nimage,obj.nCoils);
                    if obj.is2DPhaseMaps
                        P_tmp = col(repmat(obj.phaseMaps(:,:,1,indx),[1, obj.A.G.Nz]));
                    else
                        P_tmp = col(obj.phaseMaps(:,:,:,indx));
                    end
                    
                    if obj.A.G.isMasked
                        P_tmp = P_tmp(obj.A.G.mask);
                    end
                    
                    for cc = 1:obj.nCoils
                        Psense(:,cc) = obj.sen(:,cc).*exp(1j*P_tmp);
                    end

                    % Zhuoyu added T correction
                    T_indx = obj.T(:, :, indx);
                    big_number = 1e3;
                    T_indx(2,1) = T_indx(2,1) / big_number;
                    T_indx(1,2) = T_indx(1,2) / big_number;
                    T_indx(1,4) = T_indx(1,4) / big_number;
                    T_indx(2,4) = T_indx(2,4) / big_number;
                    T_indx(3,4) = T_indx(3,4) / big_number;

                    T_indx_inv = inv(T_indx);
                    tform = affine3d(T_indx_inv.');
                    disp(tform.T);
                    %disp(size(T_indx_inv));
                    x_warped_3d = imwarp(x_embed, tform, 'OutputView', imref3d([120,120,4]),'Interp','linear');
                    x_warped = x_warped_3d(obj.A.G.mask);
                    disp(size(x_warped));
                    %x_warped = x_warped_3d(:);
                    
                    
                    S = sense(obj.A_Array{indx},squeeze(Psense));
                    %Zhuoyu checking size
                    disp(size(S)); %143716*45216
                    disp(indx); %1-4
                    disp(size(x)); %45216*1
                    
                    y(:, indx) = S*x_warped;
                    
                end
                
                y = y(:);
                
            else % Adjoint operation
                
                y = zeros(obj.Nimage, obj.nShots);
                x = reshape(x, [obj.ShotLength*obj.nCoils, obj.nShots]);
                for indx = 1:obj.nShots
                    Psense = zeros(obj.Nimage, obj.nCoils);
                    
                    if obj.is2DPhaseMaps
                        P_tmp = col(repmat(obj.phaseMaps(:,:,1,indx),[1, obj.A.G.Nz]));
                    else
                        P_tmp = col(obj.phaseMaps(:,:,:,indx));
                    end
                    
                    if obj.A.G.isMasked
                        P_tmp = P_tmp(obj.A.G.mask);
                    end
                    
                    for cc = 1:obj.nCoils
                        Psense(:,cc) = obj.sen(:,cc).*exp(1j*P_tmp);
                    end
                    
                    S = sense(obj.A_Array{indx},squeeze(Psense));

                    % y(:,indx) = S'*x(:,indx);
                    %Zhuoyu added T correction
                    y_uncorr = S'*x(:,indx);
                    %y_uncorr_3d = reshape(y_uncorr,[120, 120, 4]); %reshape to 3D
                    y_uncorr_3d = embed(y_uncorr, obj.A.G.mask);

                    T_indx = obj.T(:, :, indx);
                    disp(T_indx);
                    big_number = 1e3;
                    T_indx(2,1) = T_indx(2,1) / big_number;
                    T_indx(1,2) = T_indx(1,2) / big_number;
                    T_indx(1,4) = T_indx(1,4) / big_number;
                    T_indx(2,4) = T_indx(2,4) / big_number;
                    T_indx(3,4) = T_indx(3,4) / big_number;

                    %T_indx(1,1) = T_indx(1,1) * big_number;
                    %T_indx(2,2) = T_indx(2,2) * big_number;

                    %T_indx = [1, 0, 0;
                    %          0, 1, 0;
                    %          0, 0, 1;]; %all zero 2D version
                    T_indx_inv = inv(T_indx);%inverse
                    disp(T_indx_inv);
                    %T_indx = T_indx';
                    %disp(T_indx);
                    tform = affine3d(T_indx_inv'); %conjugate transpose, cannot do tform on conjugate transpose
                    disp(tform.T);
                    y_warped_3d = imwarp(y_uncorr_3d, tform, 'OutputView', imref3d([120,120,4]),'Interp','linear');
                    %y_warped = y_warped_3d(:);
                    y_warped = y_warped_3d(obj.A.G.mask);
                    disp(size(y_warped));
                    y(:,indx) = y_warped;

                end
                y = sum(y,2);
            end
            
        end
        
    end
end

