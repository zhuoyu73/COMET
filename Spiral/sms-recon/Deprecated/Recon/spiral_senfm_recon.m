function spiral_senfm_recon(folder,recon_type,lImagetype,PreWhiten)
% INPUTS:
%   folder     - location of folder containing raw senfm file relative to
%                currect path
%   recon_type - 1: Sen map only
%                2: FM only (note this requires sen map to already have
%                already been created, except for single coil).
%                3: Sen and FM
%   lImagetype - Parameter to specify what type of image is being
%                reconstructed. Different object types can use different
%                smoothing and masking options.
%                1: full brain human
%   PreWhiten  - 1: PreWhiten Channels
%                0: Don't PreWhiten Chanels
%                TODO: add more options, PIG, single slice, leg...


% Anh: 09/12/10
% Use min-var phase error correction
% Support flexible echoes and spiral in/out readouts
% Joe:2011_01_05 make changes in order to automate recon
% Joe:2011_06_07 Convert to full brain with multislab to match sequence.
%                Only handles single echo.
% Joe:2011_10_06 Convert to use GPU, must use FM and Sense
% Joe:2011_11_10 Shift Image to center of the field of view before recon
% Joe:2011_11_10 Add recon type feature
% Joe:2011_12_06 Navigator Recon can use SENSE
% Joe:2011_12_06 Enable Multiple GPUs
% Joe:2012_06_05 Enable reconstruction of additional image sizes
% Joe:2012_06_05 Converted to use the general spiral sequence
% Joe:2012_06_18 Allow to use spiral inout
% Joe:2012_10_10 Allow for the use of averages
% Joe:2012_10_30 Modified code to be used to create sen and field maps
% Joe:2012_11_19 Modified code for sequence with combined echo times
% Joe:2012_11_26 Added type option
% Joe:2014_02_24 use fast_mr_v2 and sense svd
% Joe:2014_12_17 Added more input options
%                customized to deal with the new "spiral_senfm" sequence
% Alex:2014_01_20 Adding noise correlation matrix calculation and prepping
%                 to deal with the extra noise scan added to the
%                 sprial_senfm sequence
% Alex:2014_01_21 Cleaning up the noise correlation and pre-whiteing code
%                 to work properly whether or not you are using the
%                 head matrix channels to coils code or not
% Joe: 2015_07_24 Switching to solver_pwls_pcg for improved regularization.
%                 Modification to SENSE formulation to weight based on SOS
%                 image.
% Joe: 2015_07_29 mask addeded to sensitivity estimation for excluding
%                 voxels with large amounts of phase variation
% Joe: 2015_09_30 merging changes for noise correlation.
%                 Updated sense map estimation. Use histogram for masking.
%                 Remove channels to coils option


%suppress warnings for printing to screen
%#ok<*PRTCAL>

%% Check for irt by looking for irt folder
if (~exist('irt','dir'))
    error('It looks like the IRT is not setup! This recon requires the irt.');
end

if nargin < 4
    PreWhiten = 0;
end

%% Initialize recon
init_dir = pwd;
if ~isempty(folder)
    dir_name = strcat(pwd,'/',folder);
else
    dir_name = pwd;
end
eval(sprintf('cd %s',dir_name))
if ~(isempty(dir('*.dat')))
    fname_d = dir('*.dat');
    fname_t = fname_d(1).name;
    fname_t = fname_t(1:end-4);
else
    sprintf('Did not find a data file to reconstruct \n')
    return
end
file_name = strcat(dir_name,'/',fname_t,'.dat');
clear fname_t fname_d

if ~exist(strcat(folder,'recon_info'),'dir')
    mkdir recon_info
end

if ~exist(strcat(folder,'raw_data'),'dir')
    mkdir raw_data
end

if exist('ascconv.mat','file')
    load ascconv.mat
else
    parsasc_full((file_name));
    load ascconv.mat
end

nc = size(ascconv.asCoilSelectMeas.asList,2);
D = ascconv.sSliceArray.asSlice(1).dReadoutFOV/10;
N = ascconv.sKSpace.lBaseResolution;
nl = ascconv.sWiPMemBlock.alFree(2);
nl_design = ascconv.sWiPMemBlock.alFree(1);
nslab = ascconv.sGroupArray.asGroup.nSize;
nla_design = nl_design;
nla = nl;
nsl = 1; %this is not 3D
N_recon = N;
nsl_recon = nsl;
nImages = ascconv.lContrasts;
[read_shift, phase_shift] = image_shift( ascconv);
lpts_before_echo = 8;
% lpts_before_echo = 0;

% K-space trajectory
if exist('recon_info/kspace_info.mat','file')
    load recon_info/kspace_info.mat
    fprintf(1,'K-space information loaded\n')
else
    fprintf(1,'Setting up K-space information\n')
    gamp = 2.2; %   was always 2.4
    gslew = 140; %200;  % was 120
    gts = 10e-06;
    tsamp = 5e-06;
%     tsamp = 4e-6;
    fprintf(1,' -Calculating k-space trajectory\n')
    [Gx,Gy,kxi,kyi,~,~] = genspi(D,N,nla_design,gamp,gslew,gts);
    kxt=interp1(0:gts:gts*length(kxi)-gts,kxi,0:tsamp:gts*length(kxi)-tsamp)';
    kyt=interp1(0:gts:gts*length(kyi)-gts,kyi,0:tsamp:gts*length(kyi)-tsamp)';
    
    nk = length(kxt)-2;
    kx = zeros(nk,nla);
    ky = zeros(nk,nla);
    kxo = kxt(1:nk);
    kyo = kyt(1:nk);
    ng = length(Gx);
    gx = zeros(ng,nla);
    gy = zeros(ng,nla);
    %rotate matrix for proper orientation
    sprintf('Performing %d rotations',nla)
    phi = 2*pi/nla;
    for ii = 0:(nla-1)
        ang_rot = phi*(ii-(nla-1)*floor(ii/nla));
        kx(:,ii+1) = kxo*cos(ang_rot) + kyo*sin(ang_rot);
        ky(:,ii+1) = kyo*cos(ang_rot) - kxo*sin(ang_rot);
        gx(:,ii+1) = Gx*cos(ang_rot) + Gy*sin(ang_rot);
        gy(:,ii+1) = Gy*cos(ang_rot) - Gx*sin(ang_rot);
    end
    kx = kx(:);
    ky = ky(:);
    gx = gx(:);
    %     gy = gy(:);
    
    kx = reshape(kx,length(kx(:))/nla,nla);
    ky = reshape(ky,length(ky(:))/nla,nla);
    kz = 0;
    ww =  weight_vor(col(kx),col(ky),nla,1);    % decimated DCF
    ww(ww>2)=2;
    nro = length(kx(:))/nla;
    Nro_adc = 1024;
%     Nro_adc = 640;
    
    N_adc = ceil((nro+lpts_before_echo)/Nro_adc);
    % Two times oversampling in readout direction, real and imaginary
    lgadc = ceil(length(gx)/nla*(gts/tsamp));
    nroa = 8*ceil(lgadc/8);
    save recon_info/kspace_info kx ky kz gx ww nro gts tsamp N_adc Nro_adc lgadc nroa nl_design
    fprintf(1,' -K-space information saved\n')
end



%% Read in all data
%do not read data if data has already beeen read in

if ~exist('raw_data/raw_data_1_1.mat','file')
    %if ~exist(strcat(folder,'raw_data/raw_data_1_1.mat'),'file')
    fprintf(1,'Reading in Raw Data fron file \n')
    fid = fopen(file_name,'r');
    bytes_to_skip_hdr = fread(fid,1,'uint32');
    fseek(fid,bytes_to_skip_hdr,'bof');
    curvol_tmp = zeros(nro, nla, nsl, nslab, nc);
    
    % Deal with Noise scan
    
    NoiseScan = zeros(nc,Nro_adc,32);
    for jj = 1:32
        for ii = 1:nc
            sMDH = ice_read_mdh_va21(fid);
            if ~(sMDH.ushSamplesInScan == Nro_adc)
                sprintf('HELP: WRONG NUMBER OF SAMPLES')
                keyboard
            end
            [raw,~] = fread(fid,2*Nro_adc,'float32');
            NoiseScan(ii,:,jj) = raw(1:2:2*Nro_adc) +1j.*raw(2:2:2*Nro_adc);
        end
    end
    % Process the noise correlation matrix and calculating noise decorrelation matrix
    
    
   if ~exist('NoiseCorr.mat','file')
        NoiseCorr = calcNoiseCorr(NoiseScan);
        save('NoiseCorr.mat','NoiseCorr');
        fprintf('Noise Correlation Matrix saved \n')
    else
        load('NoiseCorr.mat');
        fprintf('Noise Correlation Matrix loaded \n')
    end
    
    
    %throw away data from discardeed acquisiitions
    for bb = 1:ascconv.sWiPMemBlock.alFree(4)
        for ii = 1:nslab
            for dd = 1:N_adc
                for cc = 1:nc
                    % Read the data from here.
                    sMDH = ice_read_mdh_va21(fid);
                    if ~(sMDH.ushSamplesInScan == Nro_adc)
                        sprintf('HELP: WRONG NUMBER OF SAMPLES')
                        keyboard
                    end
                    [~,~] = fread(fid,2*Nro_adc,'float32');
                end
            end
        end
    end
    
    for aa = 1:nImages
        fprintf(1,' -Reading data for TE %d \n',aa)
        for jj = 1:nsl
            for ll = 1:nl
                for ii = 1:nslab
                    for dd = 1:N_adc
                        for cc = 1:nc
                            % Recol(repmat(tt,1,nla))ad the data from here.
                            sMDH = ice_read_mdh_va21(fid);
                            if ~(sMDH.ushSamplesInScan == Nro_adc)
                                sprintf('HELP: WRONG NUMBER OF SAMPLES')
                                keyboard
                            end
                            [raw,~] = fread(fid,2*Nro_adc,'float32');
                            if (dd ==1)
                                if  (N_adc == 1)
                                    curvol_tmp(:,ll,jj,ii, cc) = raw(((1+2*lpts_before_echo):2:2*(lpts_before_echo+nro)))+1i*raw((2+2*lpts_before_echo):2:2*(lpts_before_echo+nro));
                                else
                                    curvol_tmp(1:(Nro_adc-lpts_before_echo),ll,jj,ii, cc) = raw(((1+2*lpts_before_echo):2:2*Nro_adc))+1i*raw((2+2*lpts_before_echo):2:2*Nro_adc);
                                end
                            elseif (dd < N_adc)&&(dd > 1)
                                curvol_tmp((((dd-1)*Nro_adc+1):(dd*Nro_adc))-lpts_before_echo,ll,jj,ii, cc) = raw(1:2:2*Nro_adc)+1i*raw(2:2:2*Nro_adc);
                            else
                                n_extra = mod(nro+lpts_before_echo,Nro_adc);
                                if n_extra == 0
                                    n_extra = Nro_adc;
                                end
                                curvol_tmp((((dd-1)*Nro_adc+1):(((dd-1)*Nro_adc+1)+n_extra-1))-lpts_before_echo,ll,jj,ii, cc) = raw(1:2:2*n_extra)+1i*raw(2:2:2*n_extra);
                            end
                        end
                    end
                end
            end
        end
        fprintf(1,sprintf(' -Saving data for TE %d\n',aa))
        
        
        if (ascconv.sSliceArray.ucMode == 4)%interleaved
            if mod(nslab,2) %odd number of slices
                exc_order = [1:2:nslab, 2:2:nslab];
            else %even number of slices
                exc_order = [2:2:nslab, 1:2:nslab];
            end
            
        elseif (ascconv.sSliceArray.ucMode == 1)%ascending
            exc_order = 1:nslab;
        elseif (ascconv.sSliceArray.ucMode == 2)%descending
            exc_order = nslab:-1:1;
        else
            exc_order = 1:nslab;
        end
        %Save Noise Correlation Data
        save raw_data/NoiseScan NoiseScan
        
        if (PreWhiten == 1)
            if ~exist('NoiseDecorr.mat','file')
                NoiseDecorr = calcNoiseDecorr(NoiseCorr);
                save('NoiseDecorr.mat','NoiseDecorr');
                fprintf('Noise Decorrelation Matrix (Prewhitening) saved \n')
            else
                load('NoiseDecorr.mat');
                fprintf('Noise Decorrelation Matrix (Prewhitening) loaded \n')
            end
        else
            NoiseDecorr = 1;
        end
        
        %save out slices seperately to limit file size
        for ii = 1:nslab
            curvol = permute(curvol_tmp(:,:,:,ii,:),[1 2 3 5 4]);
            if (PreWhiten == 1)
                for ll = 1:nl
                    for dd = 1:N_adc
                        curvol(:,ll,dd,:) = (NoiseDecorr*squeeze(curvol(:,ll,dd,:)).').';
                    end
                end
            end
            
            eval(sprintf('save raw_data/raw_data_%d_%d curvol',exc_order(ii),aa))
        end
    end
    fclose(fid);
    clear curvol sMDH bytes_to_skip_hdr curvol_tmp
    fprintf(1,' -Data read complete\n')
end

%% create sensitivity map
if (recon_type == 1)||(recon_type ==3)
    fprintf(1, ' -creating SENSE map\n')
    %need to save out coil images
    
    if ~exist('coil_imgs.mat','file')
        %if ~exist(strcat(folder,'coil_imgs.mat'),'file')
        %grid coil images
        coil_imgs = zeros(N_recon,N_recon,nslab,nc);
        for ii=1:nslab
            curvol = get_raw_data(ii,1);
%             curvol = get_raw_data(ii,2);
            for cc = 1:nc
                curvol(:,:,:,cc) = curvol(:,:,:,cc).*exp(-1i*(read_shift.*kx+phase_shift.*ky)*2*pi);
            end
            for cc = 1:nc
                coil_imgs(:,:,ii,cc) = k2image(col(kx), col(ky),col(curvol(:,:,:,cc)),ww,N_recon,2.5);
            end
        end
        save coil_imgs coil_imgs
        fprintf(1, ' -coil images saved\n')
        if nc > 1
            sos_img = sos_image(coil_imgs);
        else
            sos_img = abs(coil_imgs);
        end
        save sos_img sos_img
    else
        load coil_imgs
        load sos_img
    end
%     nc = 1;
%     nslab = 1;
    scale_factor = max(sos_img(:));
    sos_img = sos_img/scale_factor;
    coil_imgs = coil_imgs/scale_factor;
    
    [x, y] = meshgrid(1:N_recon,1:N_recon);
    mask_circ = sqrt((x-(N_recon+1)/2).^2+(y-(N_recon+1)/2).^2)<((N_recon)/2);
    mask_circ = repmat(mask_circ,[1 1 nslab]);
    N_extra = 0;
    switch lImagetype
        case 1
            %create a mask for estimate coil sensitivity
            
%         %mask objcect based on histogram of values
          [hist_counts, hist_edges]= histcounts(sos_img(:));
          ii = 1;
          while hist_counts(ii)>hist_counts(ii+1)
              ii = ii + 1;
          end
          mask = (sos_img>hist_edges(ii)).*mask_circ;
            
            %make a mask to remove data weighting from voxels with large
            %phase differences
            %TODO: rewrite with linear algebra for simplicity
            [Nx, Ny, Nz, Nc] = size(coil_imgs);
%             Rdiff = Robject(ones(size(mask(:,:,1))),'edge_type','tight','order',1,'beta',1,'type_denom','matlab','potential','quad');
            tmpMask = zeros(size(coil_imgs));
            phaseThreshhold = 0.2;
            for cc = 1:Nc
                for xx = 1:Nx
                    for yy = 1:Ny
                        for zz = 1:Nz
                            if xx >1
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) +(abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx-1,yy,zz,cc))))))<phaseThreshhold);
                            end
                            if xx < Nx
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx+1,yy,zz,cc))))))<phaseThreshhold);
                            end
                            if yy >1
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) +(abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy-1,zz,cc))))))<phaseThreshhold);
                            end
                            if yy < Ny
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy+1,zz,cc))))))<phaseThreshhold);
                            end
                            if zz >1
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) +(abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy,zz-1,cc))))))<phaseThreshhold);
                            end
                            if zz < Nz
                                tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy,zz+1,cc))))))<phaseThreshhold);
                            end
                        end
                    end
                end
            end
            if Nc == 4
                if Nz > 4
                    mask = mask.*(mean(tmpMask,4)>3.5);
                else
                    mask = mask.*(mean(tmpMask,4)>2);
                end
            elseif Nc == 32
                if Nz >4
                    mask = mask.*(mean(tmpMask,4)>2.5);
                else
                    mask = mask.*(mean(tmpMask,4)>2);
                end
            else
                mask = mask.*(mean(tmpMask,4)>2);
            end
            mask = (mask > 0);

            %smoothing terms for SENSE maps
            beta = 2e-2;
            niter = 1000; %TODO: adjust iterations based on number of voxels
            
%             init_img = 1e-2*ones(N_recon*N_recon,1);
%             Cn = C3D_sparse(ones(N_recon,N_recon))*beta;
        case 2
           %create a mask for estimate coil sensitivity
            mask = (abs(sos_img) > (0.2*max(abs(sos_img(:)))));
            
            mask = mask.*mask_circ;
            mask = (mask > 0);
            beta = 1e0;
            niter = 1000;
           
        otherwise
            %create a mask for estimate coil sensitivity
            mask = true(size(sos_img));
            Cn = 0;
            niter = 1;
            init_img = ones(N_recon*N_recon,1);
            beta = 1;
    end
    
    %mask coil images with a circular FOV
    for coilIndex= 1:nc
        coil_imgs(:,:,:,coilIndex) = coil_imgs(:,:,:,coilIndex).*mask_circ;
    end
    
    sen = zeros(N_recon,N_recon,nslab,nc);
    sen_extra = zeros(N_recon+2*N_extra,N_recon+2*N_extra,nslab,nc);

    %weight the edges to force the SENSE map to go to a certain value there
    edge_weight = mean(sos_img(:));
    data_weight = zeros(size(sen_extra(:,:,:,1)));
    data_weight((1:N)+N_extra,(1:N)+N_extra,:,1) = sos_img.*mask;
    data_weight(1,:,:) = edge_weight;
    data_weight(end,:,:) = edge_weight;
    data_weight(:,1,:) = edge_weight;
    data_weight(:,end,:) = edge_weight;

%     W = Gdiag(col(sos_img.*mask+~(mask_circ)*mean(sos_img(:))));
    if nslab < 4
        R = Robject(ones(size(sen_extra(:,:,1,1))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
        A = sensemult(ones(size(sen_extra(:,:,1,1))));
        init_img = zeros(size(sen_extra(:,:,1,1)));
        init_data = zeros(size(sen_extra(:,:,1,1)));
    else
        R = Robject(ones(size(sen_extra(:,:,:,1))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
        A = sensemult(ones(size(sen_extra(:,:,:,1))));
        W = Gdiag(data_weight(:));
        %         W = 1;
%     W = Gdiag(col(sos_img));
        init_img = zeros(size(sen_extra(:,:,:,1)));
        init_data = zeros(size(sen_extra(:,:,:,1)));
    end
    for coilIndex = 1:nc
        if nslab < 4
            for sliceIndex = 1:nslab
                init_img((1:N)+N_extra,(1:N)+N_extra) = (coil_imgs(:,:,sliceIndex,coilIndex)./sos_img(:,:,sliceIndex).*mask(:,:,sliceIndex));
                init_data((1:N)+N_extra,(1:N)+N_extra) = (coil_imgs(:,:,sliceIndex,coilIndex)./sos_img(:,:,sliceIndex).*mask(:,:,sliceIndex));
                W = Gdiag(col(data_weight(:,:,sliceIndex)));
                tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
                tmp = reshape(tmp,N+2*N_extra,N+2*N_extra);
                sen(:,:,sliceIndex,coilIndex) = tmp((1:N)+N_extra,(1:N)+N_extra,:);
            end
        else
            init_img((1:N)+N_extra,(1:N)+N_extra,:,1) = (coil_imgs(:,:,:,coilIndex)./sos_img.*mask);
            init_data((1:N)+N_extra,(1:N)+N_extra,:,1) = (coil_imgs(:,:,:,coilIndex)./sos_img.*mask);
            tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
            tmp = reshape(tmp,N+2*N_extra,N+2*N_extra,nslab);
            sen(:,:,:,coilIndex) = tmp((1:N)+N_extra,(1:N)+N_extra,:);
        end
    end
    %create an object support mask
%     mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));
     %mask objcect based on histogram of values
      [hist_counts, hist_edges]= histcounts(sos_img(:));
      ii = 1;
      while hist_counts(ii)>hist_counts(ii+1)
          ii = ii + 1;
      end
      mask = (sos_img>hist_edges(ii)).*mask_circ;
%     mask = mask.*mask_circ;
%     se1 = strel('disk',3);
%     se2 = strel('disk',5);
%     for sliceIndex = 1:nslab
%         mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se1);
%         mask(:,:,sliceIndex) = imerode(mask(:,:,sliceIndex),se2);
%         mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se2);
%     end
    mask = (mask > 0);
    
    eval(sprintf('cd %s', init_dir))
    save sen sen
    save mask mask
    eval(sprintf('cd %s',dir_name))
end

%% create field maps
if (recon_type == 2)||(recon_type ==3)
    
    fprintf(1, ' -creating field map\n')
    
    field_images = zeros(N_recon,N_recon,nslab,nImages);
    coil_rank = min(nc,3); %this data is fully sampled so we do not need a high coil rank
    if ~exist('img.mat','file')
        %if ~exist(strcat(folder,'img.mat'),'file')
        if nc > 1
            eval(sprintf('load %s/sen',init_dir))
        end
        
        eval(sprintf('load %s/mask',init_dir))
        mask = true(size(mask));
        
        niter =15;
        imginit = zeros(N_recon, N_recon, nsl_recon);
        
        for ii=1:nslab
            A = fast_mr_v2(col(kx), col(ky), 0,  D, N_recon,N_recon,1, 2*N_recon,2*N_recon,1, 5, 0, 0, 0, 0, 1,0,0,logical(mask(:,:,ii)));
            R = Robject(mask(:,:,ii),'edge_type','tight','order',1,'beta',0,'type_denom','matlab','potential','huber','delta',0.1);
            if nc > 1
                sen_tmp = reshape(sen(:,:,ii,:),N_recon*N_recon,nc);
                S = sense_svd(A ,sen_tmp(col(mask(:,:,ii)),:),coil_rank);
            end
            for aa=1:nImages
                curvol = get_raw_data(ii,aa);
                for cc = 1:nc
                    curvol(:,:,:,cc) = curvol(:,:,:,cc).*exp(-1i*(read_shift.*kx+phase_shift.*ky)*2*pi);
                end
%                 keyboard
                if nc > 1
                    c_tmp = embed(solve_pwls_pcg(col(imginit(mask(:,:,ii))), S, 1, prepData(S,col(curvol)), R, 'niter', niter),mask(:,:,ii));
                else
                    c_tmp = embed(solve_pwls_pcg(col(imginit(mask(:,:,ii))), A, 1,curvol(:), R, 'niter', niter),mask(:,:,ii));
                end
                field_images(:,:,ii,aa)  =reshape(c_tmp,N_recon,N_recon);
            end
        end
        img = field_images;
        save img img
    else
        load img
        field_images = img; %#ok<NODEF>
    end
    
    TEs = ascconv.alTE;
    
    switch lImagetype
        case 1
            % we will use a much stricter mask for field map estimate
            % this is due to bad estimates in areas of low signal that can corrupt
            % the field map
%             mask = (abs(img(:,:,:,1)) > (0.12*max(abs(col(img(:,:,:,1))))));
%             se1 = strel('disk',8);
%             se2 = strel('disk',8);
%             for sliceIndex = 1:nslab
%                 mask(:,:,sliceIndex) = imerode(mask(:,:,sliceIndex),se1);
%                 mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se2);
%             end
            mask = (abs(img(:,:,:,1)) > (0.2*max(abs(col(img(:,:,:,1))))));
            se1 = strel('disk',8);
            se2 = strel('disk',16);
            for sliceIndex = 1:nslab
                mask(:,:,sliceIndex) = imerode(mask(:,:,sliceIndex),se1);
                mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se2);
            end
            mask = (mask > 0)&(abs(img(:,:,:,1))>0);
        case 2
            mask = true(size(img(:,:,:,1)));
        otherwise
            mask = true(size(img(:,:,:,1)));
    end
    
    if nslab < 4
        FM = zeros(N_recon, N_recon, nslab);
        for sliceIndex = 1:nslab
%             keyboard
            FM(:,:,sliceIndex) = -mri_field_map_reg(squeeze(field_images(:,:,sliceIndex,1:nImages)), TEs(1:nImages)*1e-6,'l2b',-3,'mask',(mask(:,:,sliceIndex)>0));
        end
    else
        FM = -mri_field_map_reg(squeeze(field_images(:,:,:,1:nImages)), TEs(1:nImages)*1e-6,'l2b',-3,'mask',(mask>0));
    end
    
    eval(sprintf('cd %s', init_dir))
    save FM FM
    eval(sprintf('cd %s',dir_name))
    
end

%% Cleanup code
eval(sprintf('cd %s', init_dir))
end

function raw = get_raw_data(slab_num, image_num)
eval(sprintf('load raw_data/raw_data_%d_%d',slab_num,image_num));
raw = curvol;
clear curvol
end