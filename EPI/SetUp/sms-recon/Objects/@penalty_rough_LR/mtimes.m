 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
	error empty
end

v = a.v;
R = a.R;
ntp = a.ntp;
Cn = a.C;

if ~a.is.transpose 
    vi_tmp = reshape(vi,[],R);
    % uv = vi_tmp*v;
    uv = vi_tmp;
    for kk = 1:R %1:ntp
        vo_tmp(:,kk) = Cn*uv(:,kk);
    end
    % vo = vo_tmp/v;
    vo = vo_tmp(:);
else        
    vi_tmp = reshape(vi,[],R);
%     uv = vi_tmp/v;
    uv = vi_tmp;
    for kk = 1:R %1:ntp
        vo_tmp(:,kk) = Cn'*uv(:,kk);
    end
    % vo = vo_tmp*v;
    vo = vo_tmp(:);
    
end





