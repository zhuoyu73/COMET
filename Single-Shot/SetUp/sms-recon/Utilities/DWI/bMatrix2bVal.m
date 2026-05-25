function [bVals, bVecs] = bMatrix2bVal(bMatrix)
    %bMatrix2bVal Converts full b matrix to b value and b vector (q-space
    %formulation). 


    sizesB = size(bMatrix);
    if sizesB(1) ~= 3
        error 'BMatrix must be 3 x 3!'
    end
    
    if sizesB(2) ~= 3
        error 'BMatrix must be 3 x 3!'
    end
    
    for ii = 1:sizesB(3)
        [V,D] = eig(bMatrix(:,:,ii));
        
        [~, idx] = max(diag(D));
        bVals(ii) = D(idx,idx);
        bVecs(:,ii) = V(:,idx);
    end
        
end

