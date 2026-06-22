function ob = CGobj(kx, ky, kz, D, Nx, Ny, Nz, nshots, P, s_map,tt,we,L,mask)
%function ob = sense(A,sen)
%	Construct MRI object, which can do Ax and A'y operations
%  sen is nptsxncoils

% This is an attempt to make MRI object valid for sensitivity
%   encoded runs also
% Max no. of element objs 20 (previously it is 16)

%	default object


energyLevel = 0.9;

ob.n = 0;
ob.smap = 0;
ob.phase = 0;
ob.is_masked = 0;
ob.mask = [];
ob.nc = 1; %number of coils
ob.coil_rank = 1;
ob.A = {};
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
ob.A_FieldCorrection = fast_mr_v2();
ob.n = 1;
ob.VS = 1;
ob.V = 1;

if nargin == 0
	ob = class(ob, 'CGobj');
	return
end

kx = reshape(col(kx), length(col(kx))/nshots, nshots);
ky = reshape(col(ky), length(col(ky))/nshots, nshots);
kz = reshape(col(kz), length(col(kz))/nshots, nshots);

%fill object
    ob.mask = mask;
    ob.nc = size(s_map, 2);
    %ob.coil_rank = coil_rank;
    [~,sig,V] = svd(s_map,0);
    ob.V = V;
    sig = diag(sig);
    coil_rank = chooseCoilRank(sig,energyLevel);
    disp(['CGobj coil rank selection returned ' num2str(coil_rank) ' virtual coils based on ' num2str(energyLevel) ' energy level.']);
    ob.coil_rank = coil_rank;
%     if ob.nc == ob.coil_rank
%         ob.V = eye(ob.nc);
%     else
%         [~,~,V] = svd(s_map,0);
%         ob.V = V;
%     end
    VS = s_map*ob.V(:,1:coil_rank);
    ob.VS = VS(ob.mask,:);
    ob.n = nshots;
%     Fc = fast_mr_v2(kx(:,1), ky(:,1), kz(:,1),  D, Nx,Ny,Nz, 2*Nx,2*Ny,2*Nz, 5, tt, we, 0, L, 1,0,0,logical(mask)); %only prep time segmentation for first
    ob.A_FieldCorrection = fast_mr_v2(kx(:,1), ky(:,1), kz(:,1),  D, Nx,Ny,Nz, 2*Nx,2*Ny,2*Nz, 5, tt, we, 0, L, 1,0,0,logical(mask));
%     Fc1 = Fc;
%     ob.A1 = Fc;
    for idx = 1:nshots
        % the correct size of we is NxN and tt is: nk*1.
%         Fc = fast_mr_3D(kx(:,idx), ky(:,idx),kz(:,idx), D, Nx, Ny, Nz,1,1, 2*Nx, 2*Ny, 2*Nz, 5, tt, we, 0, L, 1); 
        Fc = fast_mr_v2(kx(:,idx), ky(:,idx), kz(:,idx),  D, Nx,Ny,Nz, 2*Nx,2*Ny,2*Nz, 5, tt, 0, 0, 0, 1,0,0,logical(mask));
        Fc = removemask(Fc);
        Fc = copy_sn(Fc,0);
        
        
%         Fc = copy_timeSegmentationParams(Fc,Fc1);
        
%         Fc.we = we;
%         Fc.int = int_tmp;
%         Fc.L = L;
        ob.A{idx} = Fc; 

    end
    clear Fc
    ob.smap = s_map(ob.mask,:);
%     ob.phase = P(ob.mask,:);	
    ob.phase = P;
    ob.is.empty	= logical(0);
    ob.is.transpose = logical(0);
    ob = class(ob, 'CGobj');
end

function coilRank = chooseCoilRank(sig,energyLevel)
sizeSig = size(sig);
nc = sizeSig(1);
for ii = 1:nc
    if ((sum(sig(1:ii) - sig(end))) > energyLevel*sum(sig - sig(end)))
        break;
    end
end
coilRank = ii;
end