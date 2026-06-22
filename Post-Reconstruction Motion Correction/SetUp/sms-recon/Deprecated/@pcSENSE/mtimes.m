function prod = mtimes(Cobj, v)

nobj = Cobj.n;
if (nobj == 0)
	error empty
end
% s_map = Cobj.smap;     % s_map is of size N*N x nc
sen = Cobj.VS;
nc = size(sen, 2);
P = Cobj.phase;        % P is of size N*N x nobj x nc
ns = size(Cobj.A_FieldCorrection);
nx = ns(2);	            % no. of column. (image voxels)
ny = ns(1);	            % no. of row.
mask = Cobj.mask;
%Psense = zeros(nx,nc);
bPInput2D = (numel(mask)~=size(P,1));
if (~Cobj.is.transpose)
	tmp_prod = zeros(ny*nc, nobj);   % before it is just nobj*ny
    for idx = 1:nobj
        Psense = zeros(nx,nc);

        if bPInput2D
            P_tmp = col(repmat(angle(P(:, idx)),[1 numel(mask)/size(P,1)]));
        else
        	P_tmp = angle(P(:, idx));
        end
        P_tmp = P_tmp(mask);
        for cc = 1:nc
%             Psense(:, cc) = sen(:, cc).*exp(1j*angle(P(:, idx)));
            Psense(:, cc) = sen(:, cc).*exp(1j*P_tmp);
        end
        A = copy_timeSegmentationParams(Cobj.A{idx},Cobj.A_FieldCorrection);
        A = addmask(A,mask(:));
        A = copy_sn(A,Cobj.A_FieldCorrection);
        S=sense(A,squeeze(Psense));  
%         eval(sprintf('A=sense_v2(Cobj.A%d,squeeze(Psense(Cobj.mask,:)));',idx)); 
%         tmp_prod_tmp = A*v;
%         tmp_prod(:, idx) = prepData(A,tmp_prod_tmp);
        tmp_prod(:, idx) = S*v;
    end
    prod = col(tmp_prod);
else
    prod_tmp = zeros(nx, nobj);
    vi = reshape(v,[ny*nc nobj]);
    for idx = 1:nobj
        Psense = zeros(nx,nc);

        if bPInput2D
            P_tmp = col(repmat(angle(P(:, idx)),[1 numel(mask)/size(P,1)]));
        else
        	P_tmp = angle(P(:, idx));
        end
        P_tmp = P_tmp(mask);
        for cc = 1:nc
            Psense(:, cc) = sen(:, cc).*exp(1j*P_tmp);
%             Psense(:, cc) = sen(:, cc).*exp(1j*angle(P(:, idx)));
        end  
%         eval(sprintf('A=sense_v2(Cobj.A%d,squeeze(Psense));',idx));    
        A = copy_timeSegmentationParams(Cobj.A{idx},Cobj.A_FieldCorrection);
        A = addmask(A,mask(:));
        A = copy_sn(A,Cobj.A_FieldCorrection);
        S=sense(A,squeeze(Psense));  
%         eval(sprintf('A=sense(Cobj.A%d,squeeze(Psense));',idx));
        prod_tmp(:,idx) = S'*vi(:,idx);
    end   
    prod = sum(prod_tmp,2);
end
