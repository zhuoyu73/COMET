function X = C3D_sparse(mask)
%function X = C3D_sparse(mask)
% Creates a roughness penalty matrix for the mask image.
% mask contains the penalty for each point in the image
% X is a sparse matrix containing penalties
        
    [szx szy szz] = size(mask);
    
    C_tmp = ones(szx,szy,szz);
    C_tmp2 =spdiag(col(C_tmp));
    mask1= spdiag(col(mask));
    
    C1 =C_tmp;
    C2 =C_tmp;
    C3 =C_tmp;
    C4 =C_tmp;
    C5 =C_tmp;
    C6 =C_tmp;
    
    %Set edges with no neighbors equal to zero
    C1(szx,:,:) = 0;
    C2(:,szy,:) = 0;
    C3(:,:,szz) = 0;
    C4(1,:,:) = 0;
    C5(:,1,:) = 0;
    C6(:,:,1) = 0;
    
    %shift matrices to specify a neighbor
    C1shift =circshift(C_tmp2,[-1 0]);
    C2shift =circshift(C_tmp2,[-szy 0]);
    C3shift =circshift(C_tmp2,[-szy*szx 0]);
    C4shift =circshift(C_tmp2,[1 0]);
    C5shift =circshift(C_tmp2,[szy 0]);
    C6shift =circshift(C_tmp2,[szy*szx 0]);
    
    %C-tmp2 - C4 gives a C matrix with neighbors for all points
    %multiplying by C1 removes neighbors on the edges that do not have
    %neighbors
    %mask1 contains determines the weighting for each point
    C1diff=mask1*spdiag(col(C1))*(C_tmp2-C1shift);
    C2diff=mask1*spdiag(col(C2))*(C_tmp2-C2shift);
    C3diff=mask1*spdiag(col(C3))*(C_tmp2-C3shift);
    C4diff=mask1*spdiag(col(C4))*(C4shift-C_tmp2);
    C5diff=mask1*spdiag(col(C5))*(C5shift-C_tmp2);
    C6diff=mask1*spdiag(col(C6))*(C6shift-C_tmp2);
    
     %concatanate all dimensions together
    %unique removes duplicate entries
    C_all = [C1diff; C2diff; C3diff; C4diff; C5diff; C6diff];
    %TODO: remove duplicate weighting from inside object
%     X=unique(C_all,'rows');
    X = C_all;
end

