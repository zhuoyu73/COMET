classdef NUFFT_GPU
    %NUFFT  - Convenience class for creating an NUFFT transform object
    % Replaces Brad Sutton's Gnufft and Gnufft_v2 classes.
    % Syntax:  G = NUFFT(kx,ky,kz,Nx,Ny,Nz,varargin)
    %
    % Inputs:
    %    kx      - N-D Array containing kspace in unitless [-(Nx-1)/2,Nx/2]
    %              dimensions
    %    ky      - N-D Array containing kspace in unitless [-(Ny-1)/2,Ny/2]
    %              dimensions
    %    kz      - N-D Array containing kspace in unitless [-(Nz-1)/2,Nz/2]
    %              dimensions. May be array of zeros in 2D recons.
    %    Nx      - Scalar containing size of reconstructed image in x
    %    Ny      - Scalar containing size of reconstructed image in y
    %    Nz      - Scalar containing size of reconstructed image in z, for
    %              a 2D reconstruction or operation, this must be 1.
    %   
    % Optiona Name-Value Pair Inputs:
    %    'mask'  - Image domain support mask of logicals to restrict the
    %              transform only to be computed at the true points.
    %    'VoxelBasis' - 'delta'  (Default) Compute the transform assuming
    %                            each voxel is represented in image space 
    %                            by a delta function.
    %                   'boxcar' Compute the transform by assuming each
    %                            voxel is represented in image space by a
    %                            boxcar function, requiring a sinc
    %                            weighting to be applied in data/kspace
    %                            domain
    %    'InterpMethod' - 'sparse' (Default) Computing the NUFFT using a
    %                               sparse matrix-vector multiplication 
    %                               for the spatial interpolator step. This
    %                               sparse matrix can be quite large, 
    %                               requiring very high amounts of memory 
    %                               for large problems. This method is
    %                               multithreaded for the spatial
    %                               interpolation step, which is usually
    %                               the most computationally costly step in
    %                               the computation of the NUFFT. Use
    %                               'table' if you run out of memory when
    %                               initalizing the NUFFT object.
    %                   - 'table'   This method uses a compiled mex file to
    %                               perform the spatial interpolation step.
    %                               It requires very little memory even for
    %                               very large problems, but is single
    %                               threaded which can be slow for large
    %                               problems.
    %                        
    % Outputs:
    %    G     - NUFFT class object with overloaded mtimes operator for
    %            computing the forward and adjoint NUFFT operations through
    %            A*x and A'*x.
    %
    % Example:
    %   TBD - AMC
    %
    %
    % Other m-files required: IRT nufft implementation and private directories
    %   TimeSegmentation.m, NUFFT.m
    %   
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Alex Cerjanic (1.0)
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % 4-Sep-2017; Last revision: 4-Sep-2017
    
    properties
        Kx;
        Ky;
        Kz;
        Nx;
        Ny;
        Nz;
        st;
        dims;
        isTranspose = false;
        isMasked = false;
        mask = [];
        SincWeighting = false;
        VoxelBasisWeights = 1;
        InterpMethod = '';

        % ---- PATCH START: GPU/pluggable backend awareness ----
        UseGPU = false; % try to use GPU backend if available
        Backend = 'cpu'; % 'cpu' | 'cufinufft' | 'gpunufft'
        fwdHandle = []; % function handle for forward NUFFT
        adjHandle = []; % function handle for adjoint NUFFT
        Precision = 'single'; % 'single' (recommended) or 'double'
        % ---- PATCH END ----
    end
    
    methods(Static, Access=private)
        function key = makeKey(Nx,Ny,Nz,k,prec) % use the same plan if the parameters are the same
            % trajectory hash 
            key = sprintf('Nx%d_Ny%d_Nz%d_Npts%d_sum%.5g_%s', ...
                Nx,Ny,Nz,size(k,2),sum(k(:)),prec);
        end
        function cache = getPlanCache()
            % Persistent cache, only defined once
            persistent planCache
            if isempty(planCache)
                planCache = containers.Map();
            end
            cache = planCache;
        end
        function clearPlanCache()
            cache = NUFFT_GPU.getPlanCache();
            remove(cache, keys(cache));
            %fprintf('gpuNUFFT plan cache cleared.\n');
        end
    end
    methods
        function obj = NUFFT_GPU(kx,ky,kz,Nx,Ny,Nz,varargin)
            
            J = 5;
            Jz = 5;
            K = 2;
            %% Deal with Inputs
            obj.Kx = kx;
            obj.Ky = ky;
            obj.Kz = kz;
            
            obj.Nx = double(Nx);
            obj.Ny = double(Ny);
            obj.Nz = double(Nz);
            InterpMethodString = 'sparse';
            VoxelBasis = 'delta';
            mask = [];
            p = inputParser();
            
            addOptional(p,'InterpMethod',InterpMethodString);
            addOptional(p,'mask',mask);
            addOptional(p,'VoxelBasis',VoxelBasis);
            addOptional(p,'GridOversampling',K);
            % ---- PATCH START: new options ----
            addOptional(p,'UseGPU',false);
            addOptional(p,'Backend','cpu'); % 'cpu' | 'cufinufft' | 'gpunufft'
            addOptional(p,'Precision','single');
            % ---- PATCH END ----
            parse(p,varargin{:});
            
            obj.mask = [];
            if ~isempty(p.Results.mask)
                obj.mask = logical(gather(p.Results.mask));
                obj.isMasked = true;
            else
                obj.isMasked = false;
            end
            %obj.mask = p.Results.mask;
            obj.InterpMethod = p.Results.InterpMethod;
            VoxelBasis = p.Results.VoxelBasis;
            K = p.Results.GridOversampling;
            % -----PATCH-----
            obj.UseGPU = logical(p.Results.UseGPU);
            obj.Backend = lower(p.Results.Backend);
            obj.Precision = p.Results.Precision;

            if ~isempty(obj.mask)
                obj.isMasked = true;
            end
            
            %{
            %% Prepare the NUFFT - initialize the args cell for nufft_init_v2
            if(Nz > 1)

                
                fprintf('Performing 3D NUFFT...\n');
                if strcmp(obj.InterpMethod,'table')
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'table', 2^10, 'minmax:kb'};
                    fprintf('Using table lookup...\n');
                else
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'minmax:kb', 2^10};
                    fprintf('Using sparse interpolator...\n');
                    % Estimate the size of sparse interpolator matrix
                    nnz = J*J*Jz*size(obj.Kx(:),1);
                    n = prod([K*obj.Ny,K*obj.Nx,K*obj.Nz]);
                    sparseMatrixSize = 2*(max(nnz,1)*(8+8) + (n+1) * 8);
                    sparseMatrixSizeGiB = sparseMatrixSize/(2^30);
                    fprintf('Minimum size of sparse matrix interpolator = %.1f GiB\n',sparseMatrixSizeGiB)
                    if sparseMatrixSizeGiB > get_free_mem/2
                        fprintf('Sparse Matrix is too large for system memory. Switching automatically to table method!');
                        args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'table', 2^10, 'minmax:kb'};
                        fprintf('Using table lookup...\n');
                    end
                    
                end
                
            else
                if strcmp(obj.InterpMethod,'table')
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Nx],[obj.Ny/2,obj.Nx/2]*(1), 'table', 2^10, 'minmax:kb'};
                    fprintf('Using table lookup...\n');
                else
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Nx],[obj.Ny/2,obj.Nx/2]*(1), 'minmax:kb', 2^10};
                    fprintf('Using sparse interpolator...\n');
                    % Estimate the size of sparse interpolator matrix
                    nnz = J*J*size(obj.Kx(:),1);
                    n = prod([K*obj.Ny,K*obj.Nx]);
                    sparseMatrixSize = 2*(max(nnz,1)*(8+8) + (n+1) * 8);
                    sparseMatrixSizeGiB = sparseMatrixSize/(2^30);
                    fprintf('Minimum size of sparse matrix interpolator = %.1f GiB\n',sparseMatrixSizeGiB)
                    if sparseMatrixSizeGiB > get_free_mem/2
                        fprintf('Sparse Matrix is too large for system memory. Switching automatically to table method!');
                        args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Nx],[obj.Ny/2,obj.Nx/2]*(1), 'table', 2^10, 'minmax:kb'};
                        fprintf('Using table lookup...\n');
                    end
                end
            end
            %}
            %% Setup voxel basis function
            switch VoxelBasis
                case 'delta'
                    obj.SincWeighting = false;
                    obj.VoxelBasisWeights = ones(size(obj.Kx(:)));
                case 'boxcar'
                    obj.SincWeighting = true;
                    obj.VoxelBasisWeights = sinc(obj.Kx(:)/obj.Nx).*sinc(obj.Ky(:)/obj.Ny).*sinc(obj.Kz(:)/obj.Nz);
                otherwise
                    error('Unrecognized Voxel Basis function type!');
            end
            
            %% Initialize transform
            
            %obj.st = nufft_init_v2(args{:});

            % ---- PATCH START: bind backends ----
            % default CPU handles
            obj.fwdHandle = @(x) nufft(double(x), obj.st);
            obj.adjHandle = @(x) nufft_adj(double(x), obj.st);
            ScaleFFT = 1 / sqrt(double(obj.Nx)*double(obj.Ny)*double(obj.Nz));
            
            if obj.UseGPU && strcmp(obj.Backend,'gpunufft')
            % GPU
                %fprintf('Using gpuNUFFT backend only, skipping CPU nufft_init_v2...\n');
        
                % trajectory
                if obj.Nz > 1
                    % --- 3D trajectory as 3×N ---
                    k = [single(obj.Kx(:)/obj.Nx).'; ...
                         single(obj.Ky(:)/obj.Ny).'; ...
                         single(obj.Kz(:)/obj.Nz).'];              % 3×N single
                    imageDim = double([obj.Nx, obj.Ny, obj.Nz]);  % row vector
                else
                    k = [single(obj.Kx(:)/obj.Nx).'; ...
                         single(obj.Ky(:)/obj.Ny).'];                       % 2×N single
                    imageDim = double([obj.Nx, obj.Ny]);          % row vector
                end

                assert(isa(k,'single') && isreal(k), 'gpuNUFFT: k 必须是 real single');
                %fprintf('NUFFT_GPU: k 大小 = %s\n', mat2str(size(k)));
                assert(size(k,2) > 1, 'gpuNUFFT: 轨迹点数 (M) 必须大于 1, 当前 M=%d', size(k,2));
                if obj.Nz > 1
                    assert(size(k,1)==3, 'gpuNUFFT: 3D 轨迹必须是 3×M, 当前大小=%dx%d', size(k,1), size(k,2));
                else
                    assert(size(k,1)==2, 'gpuNUFFT: 2D 轨迹必须是 2×M, 当前大小=%dx%d', size(k,1), size(k,2));
                end %ZS 0928

                w = ones(size(k,2),1,'single');   
                cache = NUFFT_GPU.getPlanCache();
                key = NUFFT_GPU.makeKey(obj.Nx,obj.Ny,obj.Nz,k,obj.Precision);
                
                if isKey(cache,key)
                    plan = cache(key);
                    %fprintf('Reusing cached plan \n');
                else
                    %fprintf('Creating new gpuNUFFT plan...\n');
                    osf = 2; wg = 4; sw = 8;
                    plan = gpuNUFFT(k, w, osf, wg, sw, imageDim, []);
                    %sens = [];
                    %atomic = true; use_textures = true; balance_workload = true;
                    %plan = gpuNUFFT(k, w, osf, wg, sw, imageDim, sens, atomic, use_textures, balance_workload);
                    cache(key) = plan;
                end

        
                if obj.Nz > 1
                    shp = [obj.Ny, obj.Nx, obj.Nz];
                else
                    shp = [obj.Ny, obj.Nx];
                end

                obj.fwdHandle = @(xg) ( ...
                    ScaleFFT * ( plan * gather(cast(xg(:), obj.Precision)) ) );
                
                obj.adjHandle = @(yg) reshape( ...
                    ScaleFFT * ( plan' * gather(cast(yg(:), obj.Precision)) ), ...
                    shp);

                obj.st = []; 
                obj.Backend = 'gpunufft';
            else
                % CPU Fessler
                if (obj.Nz > 1)
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)], ...
                            [obj.Ny, obj.Nx, obj.Nz], [J,J,Jz], [K*obj.Ny,K*obj.Nx,K*obj.Nz], ...
                            [obj.Ny/2, obj.Nx/2, obj.Nz/2], 'minmax:kb', 2^10};
                else
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)], ...
                            [obj.Ny,obj.Nx], [J,J], [K*obj.Ny,K*obj.Nx], ...
                            [obj.Ny/2,obj.Nx/2], 'minmax:kb', 2^10};
                end
        
                obj.st = nufft_init_v2(args{:});
                obj.fwdHandle = @(x) ScaleFFT * nufft(double(x), obj.st);
                obj.adjHandle = @(x) ScaleFFT * nufft_adj(double(x), obj.st);
        
                obj.Backend = 'cpu';
                %fprintf('NUFFT backend: CPU (Fessler)\n');
            end
    
            if obj.isMasked
                obj.dims = [size(obj.Kx(:),1), nnz(obj.mask(:))];
            else
                obj.dims = [size(obj.Kx(:),1), obj.Nx*obj.Ny*obj.Nz];
            end

            % ---- PATCH END ----
        end
        %% Overloaded Operators
        function y = mtimes(obj, x)
            %tic
            % Deal with forward transform first
            if ~obj.isTranspose
                
                %{
                % Apply mask to input data
                if obj.isMasked
                    x = embed(x, obj.mask);
                end
                % Reshape the data to be of the right size since it comes
                % in as a column vector
                x = reshape( x, obj.Ny,obj.Nx,obj.Nz);
                % Apply the transform using the interpolator matrix or
                % table operation
                y = nufft(double(x), obj.st);
                if (obj.SincWeighting == true)
                    y = y.*obj.VoxelBasisWeights;
                end
                %}
                % ---- Forward ----
                    if obj.isMasked
                        % mask size
                        ny = double(size(obj.mask, 1));
                        nx = double(size(obj.mask, 2));
                        nz = double(size(obj.mask, 3));  
                        if isempty(nz)
                            nz = 1; 
                        end
                    
                        % check length
                        if numel(x) ~= nnz(obj.mask)
                            error('NUFFT_GPU:MaskedInputLengthMismatch', ...
                                  'numel(x)=%d, but nnz(mask)=%d.', numel(x), nnz(obj.mask));
                        end
                    
                        if isa(x,'gpuArray')
                            tmp = complex(gpuArray.zeros(ny, nx, nz, obj.Precision));
                            m   = gpuArray(obj.mask);
                        else
                            tmp = zeros([ny nx nz], 'like', x); 
                            m   = obj.mask;
                        end
                    
                        tmp(m) = x;
                        xFull  = tmp;
                    else
                        xFull = reshape(x, double(obj.Ny), double(obj.Nx), double(obj.Nz));
                    end
                    % xFull 是 [Ny Nx (Nz)] 复数图像 或 mask embed 后的满尺寸图像
                    % 对 gpuNUFFT，我们把它展平成列向量再交给 fwdHandle
                    xFull = reshape(xFull, [], 1);
                    
                    assert(numel(xFull) > 1, 'gpuNUFFT forward: imdata 只有 1 个元素。检查 mask / reshape / 维度是否正确！');
                    %ZS 0928
                    %if strcmp(obj.Backend,'cpu') && isa(xFull,'gpuArray')
                    %    xFull = cast(xFull,'single');
                    %    if ~isa(xFull,'gpuArray'), xFull = gpuArray(xFull); end
                    %    xFull = xFull(:);  % 保证列向量
                    %    assert(numel(xFull) > 1, 'gpuNUFFT forward 输入只有 1x1，检查上游数据 reshape 是否正确');
                    %end
                    y = obj.fwdHandle(xFull);
                    if obj.SincWeighting
                        y = y .* cast(obj.VoxelBasisWeights, 'like', y);
                    end
                
            else % Deal with the adjoint case
                %{
                % Apply the adjoint transform
                if (obj.SincWeighting == true)
                    x = x.*obj.VoxelBasisWeights;
                end
                y = nufft_adj(double(x), obj.st);
                if obj.isMasked
                    % Return only masked data
                    y = y(obj.mask);
                end
                %}
                % ---- Adjoint ----
                
                xK = x(:);  % kspace data 展平成列向量
                %fprintf('adjoint 输入大小 = %s\n', mat2str(size(xK)));

                % 断言：样本数 > 1
                assert(numel(xK) > 1, 'gpuNUFFT adjoint: data 只有 1 个元素。检查 data_US / dataMask / 维度是否正确！');

                %xK = x; %0928 ZS
                if obj.SincWeighting
                    xK = xK .* cast(obj.VoxelBasisWeights, 'like', xK);
                end
                if strcmp(obj.Backend,'cpu') && isa(xK,'gpuArray')
                        xK = cast(xK,'single');
                end
                %fprintf('adjoint 输入大小 = %s\n', mat2str(size(xK)));
                yImg = obj.adjHandle(xK);
            
                if obj.isMasked
                    if isa(yImg,'gpuArray')
                        m = gpuArray(obj.mask);
                    else
                        m = obj.mask;
                    end
                    y = yImg(m);
                else
                    y = yImg;
                end


            end
            % Ensure we return a column vector
            y = y(:);
            %toc
        end
        
        function obj = ctranspose(obj)
            % Handle the transpose operator for adjoint operator
            obj.isTranspose = ~obj.isTranspose;
            obj.dims = fliplr(obj.dims);
        end
        
        function y = size(obj)
            % Return the size of the object
            y = obj.dims;
        end
        
        function y2 = make2xN_realimag(yg, prec)
            yg=gather(yg(:));
            y2=cast([real(yg).'; imag(yg).'],prec);
        end
                
    end
    
end

