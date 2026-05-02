function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
    error empty
end

i = sqrt(-1);

L = size(a.int,1);

if (L == 1)
    tau = 0;
else
    tau = (max(a.tt)-min(a.tt)+eps)/(L-1);
end
TE = min(a.tt);
sizeG = size(a.G);


if ~a.is.transpose
    vi = exp(-i*a.we(:)*TE).*vi(:);
    vo_tmp = zeros(sizeG(1),L);
    for ll = 1:L
        Wo = exp(-i*col(a.we)*((ll-1)*tau));
        vo_tmp(:,ll) = a.G*(Wo.*col(vi));
    end
    aa = repmat(permute(a.int, [2 1]), [a.ttsegments, 1]);
    vo = sum(vo_tmp.*aa, 2);
    if a.flgswt
        vo = vo.*a.swt;
    end
    
else
    if a.flgswt
        vi = vi.*a.swt;  % Transpose of real sinc-weighting
    end
    vo_tmp = zeros(sizeG(2), L);
    
    for ll = 1:L
        Wo = exp(i*conj(a.we)*((ll-1)*tau));
        aa = repmat(a.int(ll,:)',[a.ttsegments 1]);
        vo_tmp(:,ll) = Wo(:).*(a.G'*(aa.*col(vi)));
    end
    vo = sum(vo_tmp, 2);
    vo = exp(i*conj(a.we(:))*TE).*vo;
    %vo = vo(:);
end





